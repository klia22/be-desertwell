# Summary of the Algorithmic Process

## The Brute‑Force Baseline

Minecraft Bedrock Edition places desert wells via a pseudo‑random lottery. The process for a given world seed (64‑bit signed) and chunk coordinates (X,Z) is:

1. Derive structure multipliers `xMul` and `zMul` from the world seed via the MT19937_32.
2. Compute a region seed  
   R = S XOR (xMul * X + zMul * Z)   (mod 2^32, signed)  
   where S is the lower 32 bits of the world seed.
3. Use R as an RNG seed to draw a PRNG MT call: `nextInt<500>() == 0` (probability 1/500).
4. If the PRNG MT call succeeds, two more PRNG calls choose an offset inside the 16×16 block, yielding final world coordinates.

We wanted to find seeds that produce a well in a **2×2 chunk section** (NW, NE, SW, SE) with specific offset ranges (0‑4 for low, 11‑15 for high), and use specific pairs of low and high to form each corner. A naive brute‑force approach would enumerate all world seeds and, for each, examine many chunk coordinates. Since the structure multipliers depend only on the lower 32 bits (structural potential), the core search space is 2^32 lower‑32 seeds. If one were to test all possible chunk coordinates in a 2^24 × 2^24 chunk grid, the total number of seed‑chunk‑corner evaluations would be 2^24 × 2^24 × 2^32 = 2^80 structure calls (which includes its own 624 MT states per each check). Tackling the problem this way, even with many optimizations and an extremely strong CPU, is astronomically impossible.

## Precomputing the Base Seeds

We observed that once the region seed R is computed from the world seed and chunk coordinates, all subsequent structure placement decisions – the 1/500 chance and the two offset draws – depend **only on R**. Since there are only 2^32 values of R, we could replace the per‑seed structure RNG testing with a **one‑time exhaustive precomputation**:

- We enumerated all 2^32 possible region seeds.
- For each, we simulated exactly what the game does: the initial structure call (`nextInt<500>() == 0`) and, if successful, the two offset draws.
- If the structure call succeeded **and** the offsets fell into the required ranges for a corner, we stored the region seed in that corner’s set.

This **base‑scanning phase** runs only once, producing four lists L_NW, L_NE, L_SW, L_SE of “good” region seeds. The size of each list is far smaller than 2^32 because it requires both the 1/500 chance and the specific offset windows. After this step, the expensive RNG logic is never needed again. The reduction is enormous: we traded the infeasible 2^80 per‑seed chunk evaluations for a single 2^32 scan (which is extremely fast, taking just 30 minutes on typical CPU), and the remaining search space becomes a 2^64 candidate check.

## Algorithm After Precomputation

With the precomputed corner sets, the problem becomes: for each world seed S (lower 32 bits), find a value U = xMul * X + zMul * Z such that four simultaneous inclusions hold:

    S XOR U               ∈ L_NW
    S XOR (U + xMul)      ∈ L_NE
    S XOR (U + zMul)      ∈ L_SW
    S XOR (U + xMul + zMul) ∈ L_SE

If we naively scanned all 2^32 possible U for each of the 2^32 seeds, the total work would be 2^64 checks. Considering rarity checks, this can be simplified to ~2^44. This is already a huge improvement from 2^80, but still far too large for a normal CPU, although extremely fast CPUs may still be able to handle it in reasonable amount of time. But there exists further optimizations.

## Adding MITM

To make the search feasible, we employ a **Meet‑in‑the‑Middle (MITM)** decomposition on the **low 22 bits** of the region seeds. We only take the lower 22 bits since just taking 22 bits approximately seperates all regions seeds into their own hash cell without taking extra memory. We do MITM and split 22 bits into:

- HIGH_BITS = 12 (4096 possible values)
- LOW_BITS = 10 (1024 possible values)

(These constants chosen are proven to be optimal for efficiency)

We indexed all precomputed region seeds by these low 22 bits into a **2D bit‑mask table** for each corner. For each corner, we constructed **4096 rows**, each row being a **1024‑bit mask** (16 × 64‑bit words) indicating which low‑10 values appear for that row.

Now, for a fixed seed S, we iterate over the **1024 possible low‑10 values of U**, denoted k. For each k, we want to find the high‑12 bits h that satisfy all four corner constraints simultaneously. This transforms the problem from a giant 4‑way membership test over all U into a series of **bit‑mask intersections** over the 4096‑element high space, one per low‑10 value.

## Carry Handling with Shift and Masks

The additive terms xMul, zMul, xMul+zMul cause carries from the low‑10 part into the high‑12 part. For a given corner, say SW, the condition is:

    S XOR (U + xMul) ∈ L_SW

If we write U_low = k and U_high = h, then the low‑10 sum k + xMulLow may overflow, adding 1 to the high‑12 part. We handled this by **precomputing two carry masks** for every possible shift amount (0…1023): a **carry‑0 mask** that marks the valid column range when no overflow occurs, and a **carry‑1 mask** for the overflow case.

For a given k, we computed a **4096‑bit mask** of valid high‑12 values for SW by:

1. Rotating the SW row corresponding to the target column `(k + xLow10) mod 1024` by the word‑shift amount (using SIMD unaligned loads and bit shifts).
2. ANDing with the carry‑0 mask.
3. If the low addition did overflow, also loading a second row (with the incremented high part) and ANDing with the carry‑1 mask, then OR‑ing the two contributions.

We did the same for NE and NW, yielding three 4096‑bit masks.

## SE Column and Final Intersection

The SE condition is simpler: it does not involve an additive shift from U; it only requires that the SE row corresponding to h XOR (…) has a set bit in column k. This is a direct column lookup: we built, for each low‑10 value, a **4096‑bit mask** indicating which high‑12 rows have a set bit in that column. (This column mask was constructed on the fly from the pre‑permuted SE rows for the current seed’s low‑10 value.)

The final candidate high‑12 values h for a given k are exactly those bits set in the intersection:

    SE_col_mask  AND  SW_mask  AND  NE_mask  AND  NW_mask

Only the bits that survive this intersection need further scrutiny.

## SIMD Acceleration

All the 4096‑bit mask operations are implemented with **AVX2 SIMD intrinsics**. A 4096‑bit mask is stored as **16 × 256‑bit registers** (`__m256i`). The key SIMD operations include:

- **Loading unaligned data** (`_mm256_loadu_si256`) to read rows from memory, and when a column shift is required, loading two adjacent 256‑bit chunks and combining them with `_mm256_srl_epi64` and `_mm256_sll_epi64` to create the shifted mask.
- **Bitwise AND / OR** across entire masks using `_mm256_and_si256` / `_mm256_or_si256`.
- **Zero‑testing** with `_mm256_testz_si256` to quickly skip rows that have no surviving candidates, avoiding the cost of bit extraction.

Additionally, we implemented a **bit‑permutation function** (`xorPermute1024AVX2`) that reorders the bits within a row according to the low‑10 value of the seed. This permutation uses a series of masked swaps and full‑register permutes (`_mm256_permute4x64_epi64`) to rearrange bits without scalar loops. We pre‑permuted the SE rows once per low‑10 value, reusing them for all seeds sharing that value, further reducing per‑seed work.

These SIMD techniques allowed us to process 256 bits per instruction, making the mask intersection loop extremely fast – the inner 4096‑row scan for a single low‑10 value completes in a few hundred CPU cycles.

## Bit Extraction and Full 32‑Bit Validation

When the intersection mask is non‑zero, we store the 4096‑bit result to a local buffer and use the scalar `__builtin_ctzll` intrinsic to extract set bits efficiently. Each set bit gives a candidate high‑12 value h. Together with k, we form the full low‑22 bits of U and compute the full 32‑bit SE base:

    baseSE = S XOR (h << 10 | k)

Because the MITM filter only uses the low 22 bits, false positives are possible. Therefore we validate every candidate against the **full 32‑bit corner sets**:

- For SW, NE, NW we use **dense bit‑sets** (256 MB arrays indexed by the 32‑bit value) – a single bit test confirms membership in O(1).
- For SE we use a **compact bucketed index** (`seBucketStart`/`seBucketCount`) that stores the full SE values grouped by their low‑22 bits, enabling fast enumeration without a second giant bit‑set.

Only if all four full‑bit checks pass do we count a candidate as a valid solution.

## Finding the Nearest Chunk Coordinates

Once a valid seed is confirmed, we still want the actual chunk coordinates (X,Z) that produce the well. This requires solving the linear Diophantine equation:

    xMul * X + zMul * Z ≡ C  (mod 2^32, signed)

where C = baseSE. We solve it using the extended Euclidean algorithm to find a particular solution. However, many integer solutions exist, and we need the one **closest to the origin**. In Minecraft Bedrock Edition, gameplay physics (player position, collision, movement) uses single‑precision floating point. At large coordinate values (e.g., beyond 2^24 blocks), the precision points extremely large (e.g, slicing at 0.5 for 2^23, 1.0 for 2^24...), causing severe jitter, falling through the world, and loss of control. Thus, only wells nearer the origin are practically visitable. We therefore minimize the **Chebyshev distance** max(|X|,|Z|) via a binary search over the integer solution space, yielding the nearest usable coordinates.

## The Full 64‑Bit Seed Extension

This entire search only scans the lower 32 bits of the world seed, because the structure multipliers (and hence well placements) depend exclusively on those bits. The upper 32 bits determine other world properties (mostly only noise properties), most notably **biomes**. To obtain a complete seed that also places a desert biome at the well location, one must later filter the candidate potentials with a biome check. This can be done efficiently using the external **cubiomes** library, which provides exact biome queries for full 64‑bit seeds. The lower‑32‑bit scan thus supplies a shortlist of potentials, which are then verified for desert biome compatibility. Still, after all of this, some in-game conditions may differ like desert wells overriding each other or caves generating overriding them.

## Compilation

We now compiled the final program with the following flags:

`-Ofast -march=native -mtune=native -flto -DNDEBUG -pthread`

These enable aggressive optimizations, target the specific host CPU, employ link‑time optimization, disable assertions, and support multithreading.

## Summary

We transformed an infeasible 2^80 brute‑force search into a practical ~2^38 scan by:

- Precomputing the region seeds over all 2^32 region seeds, removing the RNG overhead.
- Using a 12‑high / 10‑low MITM split to reduce per‑seed work from scanning all U to 1024 fast mask intersections.
- Handling carries with precomputed shift masks.
- Accelerating all 4096‑bit mask operations with AVX2 SIMD intrinsics.
- Applying final entire 32‑bit bit‑set checks to ensure zero false positives.
- Solving for the nearest chunk coordinates via Diophantine approximation, respecting Minecraft’s floating‑point limitations.

The result is a system that can scan the entire 2^32 lower‑seed space on a desktop machine (around 5 days on a standard machine), making a previously impossible task routine.
