// These functions are direct copies from ezseed.net's internal code worker-7ASCFYwv.js
class X {
    constructor(t) {
        this.rng = new $r(mo(t))
    }
    setSeed(t) {
        this.rng._seed(mo(t))
    }
    nextInt(t) {
        return t == null ? this.rng.random_int() >>> 1 : this.rng.random_int() % t
    }
    nextFloat() {
        return Math.fround(this.rng.random_int() / 4294967296)
    }
    nextDouble() {
        return this.rng.random_int() / 4294967296
    }
    nextBoolean() {
        return (this.rng.random_int() >>> 27 & 1) === 1
    }
}
class $r {
    constructor(t) {
        this.mt = new Uint32Array(624),
        this.index = 624,
        this._seed(t >>> 0)
    }
    _seed(t) {
        this.mt[0] = t >>> 0;
        for (let n = 1; n < 624; n++) {
            const o = this.mt[n - 1];
            this.mt[n] = Math.imul(1812433253, (o ^ o >>> 30) >>> 0) + n >>> 0
        }
        this.index = 624
    }
    _twist() {
        for (let i = 0; i < 624; i++) {
            const r = (this.mt[i] & 2147483648 | this.mt[(i + 1) % 624] & 2147483647) >>> 0;
            let a = (this.mt[(i + 397) % 624] ^ r >>> 1) >>> 0;
            r & 1 && (a = (a ^ 2567483615) >>> 0),
            this.mt[i] = a
        }
        this.index = 0
    }
    random_int() {
        this.index >= 624 && this._twist();
        let t = this.mt[this.index++];
        return t = (t ^ t >>> 11) >>> 0,
        t = (t ^ t << 7 & 2636928640) >>> 0,
        t = (t ^ t << 15 & 4022730752) >>> 0,
        t = (t ^ t >>> 18) >>> 0,
        t
    }
}
function mo(e) {
    return typeof e == "bigint" ? Number(BigInt.asIntN(32, e)) >>> 0 : e >>> 0
}
function ia(e, t, n) {
    const o = (e.seedLow32 ^ (Math.imul(t, e.xMul) + Math.imul(n, e.zMul) | 0)) >>> 0;
    return (o ^ e.featureKey + (o << 6 >>> 0) + (o >>> 2) - ea >>> 0) >>> 0
}
function na(e) {
    let t = -2078137563;
    for (let n = 0; n < e.length; n++)
        t = Math.imul(t, 435) ^ e.charCodeAt(n) | 0;
    return t >>> 0
}
function oa(e, t) {
    const n = typeof e == "bigint" ? e : BigInt(e)
      , o = Number(BigInt.asIntN(32, n)) | 0
      , i = new X(o);
    return {
        seedLow32: o,
        xMul: (i.nextInt() | 1) >>> 0,
        zMul: (i.nextInt() | 1) >>> 0,
        featureKey: na(t)
    }
}
function O(e, t) {
    let n = e / t | 0;
    return e % t !== 0 && (e ^ t) < 0 && n--,
    n
}
const ea = 1640531527
  , Zo = Object.freeze({
    Geode: "minecraft:overworld_amethyst_geode_feature",
    DesertWell: "minecraft:desert_after_surface_desert_well_feature",
    Fossil: "minecraft:desert_or_swamp_after_surface_fossil_feature",
    FossilDeepslate: "minecraft:desert_or_swamp_after_surface_fossil_deepslate_feature"
});
const ra = Zo.DesertWell;
function aa(e, t, n, o, i, r=null) {
    const a = oa(e, ra)
      , s = []
      , c = O(t, 16)
      , u = O(n, 16)
      , l = O(t + o - 1, 16)
      , f = O(n + i - 1, 16);
    for (let v = u; v <= f; v++)
        for (let d = c; d <= l; d++) {
            if (r && !r(d * 16 + 8, v * 16 + 8))
                continue;
            const h = new X(ia(a, d, v));
            if (h.nextInt(500) >= 1)
                continue;
            const _ = h.nextInt(16)
              , m = h.nextInt(16)
              , g = d * 16 + m
              , p = v * 16 + _;
            g < t || g >= t + o || p < n || p >= n + i || r && !r(g, p) || s.push({
                x: g,
                z: p
            })
        }
    return s
}
// Parameters: e: seed, t: startX, n: startZ, o: width, i: height
// Self reminder to run a 32 upper 32 search for biome matching
async function main() {
    console.log(aa(9288674231451659n, -1024, -1024, 2048, 2048)); // Note, the n is important and cannot be removed
}
main();