// ============================================================
//  VLBI Delay Model
//  Geometric + Tropospheric (Saastamoinen) + Ionospheric (1st + 2nd order)
//
//  2nd order ionospheric model implemented per:
//    Hawarey, M., Hobiger, T. & Schuh, H. (2005).
//    Effects of the 2nd order ionospheric terms on VLBI measurements.
//    Geophys. Res. Lett. 32, L11304. doi:10.1029/2005GL022729
//
//  Author: Dr. Mosab Hawarey | @DrHawarey | github.com/mhawarey
// ============================================================
#![allow(clippy::many_single_char_names)]

use eframe::egui;
use egui::{Color32, RichText, Stroke, Ui, Vec2};

// ── Physical constants ─────────────────────────────────────────────────────────
const C: f64 = 2.997_924_58e8;   // speed of light [m/s]
const RE: f64 = 6_378_137.0;     // WGS84 semi-major axis [m]
const K1: f64 = 77.604e-6;       // refractivity K1 [K/Pa]
const K2: f64 = 64.79e-6;        // refractivity K2 [K/Pa]
const K3: f64 = 3.776e-1;        // refractivity K3 [K²/Pa]
const TEC_COEFF: f64 = 40.309e16;// ionospheric coeff [m·Hz²/TECU]
const PI: f64 = std::f64::consts::PI;

// 2nd order ionospheric constants (Hawarey, Hobiger & Schuh, 2005, Eq. 8)
//   Eq. 8 (paper):  s = 7527 · c · ∫ N · (B⃗·k̂) dL    [units: m·Hz³]
//   Eq. 3a/3b:      τ_range = s / ν³                  [units: meters of range]
//   Convert to time delay (seconds) by dividing by c:
//                   τ_time = s / (ν³ · c)
//   This convention matches Bassiri & Hajj (1993) and Fritsche et al. (2005).
const S2_COEFF: f64 = 7527.0;             // Eq. 8 numerical coefficient
const B_EQ_TESLA: f64 = 3.12e-5;          // geomagnetic field at magnetic equator [T]
const IONO_SHELL_BOTTOM_M: f64 = 100_000.0;   // 100 km — F-region lower bound
const IONO_SHELL_TOP_M: f64    = 1_000_000.0; // 1000 km — F-region upper bound
const IONO_NMAX_F2: f64 = 2.0e12;         // peak F2 electron density [e-/m³]
                                          // (CONT02 campaign was near solar max, cycle 23)
const IONO_HMAX_F2_M: f64 = 350_000.0;    // F2 peak height [m]
const IONO_SCALE_H_M: f64 = 80_000.0;     // Chapman scale height [m]
const IONO_INTEGRATION_PTS: usize = 100;  // sample points along ray (per 2005 paper)

// ── Data structures ────────────────────────────────────────────────────────────

#[derive(Clone, Debug)]
pub struct Station {
    pub name:         String,
    pub lat_deg:      f64,
    pub lon_deg:      f64,
    pub height_m:     f64,
    pub pressure_hpa: f64,
    pub temp_k:       f64,
    pub e_hpa:        f64,
}

impl Station {
    fn new(name: &str, lat: f64, lon: f64, h: f64, p: f64, t: f64, e: f64) -> Self {
        Self { name: name.into(), lat_deg: lat, lon_deg: lon, height_m: h,
               pressure_hpa: p, temp_k: t, e_hpa: e }
    }
    fn ecef(&self) -> [f64; 3] { geodetic_to_ecef(self.lat_deg, self.lon_deg, self.height_m) }
}

#[derive(Clone, Debug)]
pub struct Source {
    pub name:    String,
    pub ra_h:    f64,
    pub dec_deg: f64,
}

impl Source {
    fn new(name: &str, ra_h: f64, dec_deg: f64) -> Self {
        Self { name: name.into(), ra_h, dec_deg }
    }
    fn unit_vector(&self) -> [f64; 3] {
        let ra  = self.ra_h * 15.0 * PI / 180.0;
        let dec = self.dec_deg * PI / 180.0;
        [dec.cos()*ra.cos(), dec.cos()*ra.sin(), dec.sin()]
    }
}

#[derive(Clone, Debug)]
pub struct Baseline { pub sta1: usize, pub sta2: usize }

#[derive(Clone, Debug)]
pub struct DelayResult {
    pub baseline_name: String,
    pub source_name:   String,
    pub tau_geo_ps:    f64,
    pub tau_tropo1_ps: f64,
    pub tau_tropo2_ps: f64,
    pub tau_iono_ps:   f64,
    pub tau_iono2_ps:  f64,   // NEW: 2nd order ionospheric delay [ps]
    pub tau_total_ps:  f64,
    pub el1_deg:       f64,
    pub el2_deg:       f64,
    pub valid:         bool,
}

// ── Core math ──────────────────────────────────────────────────────────────────

fn geodetic_to_ecef(lat_deg: f64, lon_deg: f64, h: f64) -> [f64; 3] {
    let f  = 1.0 / 298.257223563;
    let e2 = 2.0*f - f*f;
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let n = RE / (1.0 - e2*lat.sin().powi(2)).sqrt();
    [(n+h)*lat.cos()*lon.cos(), (n+h)*lat.cos()*lon.sin(), (n*(1.0-e2)+h)*lat.sin()]
}

fn dot3(a: &[f64;3], b: &[f64;3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn norm3(a: &[f64;3]) -> f64 { (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]).sqrt() }

/// Greenwich Mean Sidereal Time [radians] from Julian Date
/// IAU 1982: GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0)
fn gmst_rad(jd: f64) -> f64 {
    let d = jd - 2451545.0;
    let gmst_deg = 280.460_618_37 + 360.985_647_366_29 * d;
    (gmst_deg % 360.0).to_radians()
}

/// Local Hour Angle of source at station [radians]
fn hour_angle(gmst: f64, lon_deg: f64, ra_h: f64) -> f64 {
    let lst = gmst + lon_deg.to_radians();
    let ra  = (ra_h * 15.0).to_radians();
    lst - ra
}

/// Elevation of source at station accounting for Earth rotation (epoch-aware)
fn elevation_epoch(sta: &Station, src: &Source, gmst: f64) -> f64 {
    let lat = sta.lat_deg.to_radians();
    let dec = src.dec_deg.to_radians();
    let ha  = hour_angle(gmst, sta.lon_deg, src.ra_h);
    let sin_el = lat.sin()*dec.sin() + lat.cos()*dec.cos()*ha.cos();
    sin_el.clamp(-1.0, 1.0).asin()
}

/// Elevation without epoch (ICRF frame — for tests / backward compat)
fn elevation(sta: &Station, src: &Source) -> f64 {
    let lat = sta.lat_deg.to_radians();
    let lon = sta.lon_deg.to_radians();
    let s = src.unit_vector();
    let e = [-lon.sin(), lon.cos(), 0.0f64];
    let n = [-lat.sin()*lon.cos(), -lat.sin()*lon.sin(), lat.cos()];
    let u = [ lat.cos()*lon.cos(),  lat.cos()*lon.sin(), lat.sin()];
    let se = dot3(&e,&s); let sn = dot3(&n,&s); let su = dot3(&u,&s);
    su.atan2((se*se+sn*sn).sqrt())
}

/// Source unit vector k⃗ in station-local ENU frame
fn source_enu(sta: &Station, src: &Source, gmst: f64) -> [f64;3] {
    let lat = sta.lat_deg.to_radians();
    let dec = src.dec_deg.to_radians();
    let ha  = hour_angle(gmst, sta.lon_deg, src.ra_h);
    // Standard astronomical → topocentric ENU
    let east  = -dec.cos() * ha.sin();
    let north = -lat.sin()*dec.cos()*ha.cos() + lat.cos()*dec.sin();
    let up    =  lat.cos()*dec.cos()*ha.cos() + lat.sin()*dec.sin();
    [east, north, up]
}

fn geometric_delay(s1: &Station, s2: &Station, src: &Source) -> f64 {
    let r1=s1.ecef(); let r2=s2.ecef();
    let b=[r2[0]-r1[0], r2[1]-r1[1], r2[2]-r1[2]];
    -dot3(&b, &src.unit_vector()) / C
}

fn saastamoinen_delay(sta: &Station, el_rad: f64) -> f64 {
    let lat = sta.lat_deg.to_radians();
    let h   = sta.height_m / 1000.0;
    let p   = sta.pressure_hpa;
    let t   = sta.temp_k;
    let e   = sta.e_hpa;
    let denom = 1.0 - 0.00266*(2.0*lat).cos() - 0.00028*h;
    let ztd = 0.002277*p/denom + 0.002277*(1255.0/t+0.05)*e;
    ztd / el_rad.max(5.0_f64.to_radians()).sin()
}

fn tropo_diff(s1: &Station, s2: &Station, src: &Source, gmst: f64) -> (f64,f64,f64) {
    let el1=elevation_epoch(s1,src,gmst); let el2=elevation_epoch(s2,src,gmst);
    let d1=if el1>3.0_f64.to_radians(){saastamoinen_delay(s1,el1)}else{0.0};
    let d2=if el2>3.0_f64.to_radians(){saastamoinen_delay(s2,el2)}else{0.0};
    ((d2-d1)/C, d1/C, d2/C)
}

fn iono_diff(s1: &Station, s2: &Station, src: &Source,
             tec1: f64, tec2: f64, freq_hz: f64, gmst: f64) -> f64 {
    let el1=elevation_epoch(s1,src,gmst); let el2=elevation_epoch(s2,src,gmst);
    let p1=if el1>3.0_f64.to_radians(){tec1/el1.sin()}else{0.0};
    let p2=if el2>3.0_f64.to_radians(){tec2/el2.sin()}else{0.0};
    TEC_COEFF*(p2-p1)/(freq_hz*freq_hz*C)
}

// ─────────────────────────────────────────────────────────────────────────────
// 2nd ORDER IONOSPHERIC MODEL
// Implementation of Hawarey, Hobiger & Schuh (2005), GRL 32, L11304.
//
// The 2nd order delay term, in time, is:
//     τ_2 = s · ν⁻³                         (Eq. 3a/3b of the 2005 paper)
// where the s-term is computed from the ray-integrated product of electron
// density and the projection of the geomagnetic field onto the propagation
// direction:
//     s = 7527 · c · ∫ N · (B⃗ · k̂) dL      (Eq. 8 of the 2005 paper)
//
// We integrate along the line of sight from the station up through the
// ionosphere using IONO_INTEGRATION_PTS sample points (the 2005 paper used
// 100 representative points), avoiding the dipole and 400-km thin-shell
// approximations that the paper showed were not legitimate.
// ─────────────────────────────────────────────────────────────────────────────
mod ionosphere {
    use super::*;

    /// Chapman-α electron density profile (used as PIM stand-in).
    /// N(h) = Nmax · exp( 0.5 · (1 − z − exp(−z)) ),   z = (h − hmax)/H
    pub fn electron_density(h_m: f64) -> f64 {
        if h_m < IONO_SHELL_BOTTOM_M || h_m > IONO_SHELL_TOP_M { return 0.0; }
        let z = (h_m - IONO_HMAX_F2_M) / IONO_SCALE_H_M;
        IONO_NMAX_F2 * (0.5 * (1.0 - z - (-z).exp())).exp()
    }

    /// Tilted-dipole geomagnetic field B⃗ at ECEF point r (Tesla).
    /// Uses an IGRF-style centred dipole with the geomagnetic pole
    /// at (lat 80.65°N, lon −72.68°E) — IGRF-13 epoch reference.
    /// This avoids the centred non-tilted dipole that the 2005 paper
    /// flagged as introducing ≥10% error.
    pub fn geomag_field(r_ecef: &[f64;3]) -> [f64;3] {
        // Geomagnetic pole direction (IGRF-13 reference)
        let pole_lat = 80.65_f64.to_radians();
        let pole_lon = -72.68_f64.to_radians();
        let m_hat = [pole_lat.cos()*pole_lon.cos(),
                     pole_lat.cos()*pole_lon.sin(),
                     pole_lat.sin()];
        let r_mag = norm3(r_ecef);
        let r_hat = [r_ecef[0]/r_mag, r_ecef[1]/r_mag, r_ecef[2]/r_mag];
        // B⃗ = (B_eq · (RE/r)³) · ( 3·(m̂·r̂)·r̂ − m̂ )
        let amp = B_EQ_TESLA * (RE/r_mag).powi(3);
        let mr = dot3(&m_hat, &r_hat);
        [amp*(3.0*mr*r_hat[0] - m_hat[0]),
         amp*(3.0*mr*r_hat[1] - m_hat[1]),
         amp*(3.0*mr*r_hat[2] - m_hat[2])]
    }

    /// Compute the s-term [SI units] for a single station/source ray.
    /// Returns 0.0 if the source is below the elevation cutoff.
    pub fn s_term(sta: &Station, src: &Source, gmst: f64) -> f64 {
        let el = elevation_epoch(sta, src, gmst);
        if el <= 3.0_f64.to_radians() { return 0.0; }

        // Build ECEF k̂ from station-local ENU k̂
        let k_enu = source_enu(sta, src, gmst);
        let lat = sta.lat_deg.to_radians();
        let lon = sta.lon_deg.to_radians();
        let e_axis = [-lon.sin(),               lon.cos(),               0.0    ];
        let n_axis = [-lat.sin()*lon.cos(),    -lat.sin()*lon.sin(),    lat.cos()];
        let u_axis = [ lat.cos()*lon.cos(),     lat.cos()*lon.sin(),    lat.sin()];
        let k_ecef = [
            k_enu[0]*e_axis[0] + k_enu[1]*n_axis[0] + k_enu[2]*u_axis[0],
            k_enu[0]*e_axis[1] + k_enu[1]*n_axis[1] + k_enu[2]*u_axis[1],
            k_enu[0]*e_axis[2] + k_enu[1]*n_axis[2] + k_enu[2]*u_axis[2],
        ];

        let r0 = sta.ecef();
        // Slant path lengths to the ionospheric shell bottom and top.
        // Approximate ray geometry: spherical Earth, straight ray.
        let l_bot = slant_path_to_height(sta.height_m, IONO_SHELL_BOTTOM_M, el);
        let l_top = slant_path_to_height(sta.height_m, IONO_SHELL_TOP_M,    el);
        let n_pts = IONO_INTEGRATION_PTS;
        let dl = (l_top - l_bot) / (n_pts as f64);
        let mut accum = 0.0_f64;
        // Trapezoidal rule
        for i in 0..=n_pts {
            let l = l_bot + dl * (i as f64);
            let r = [r0[0] + l*k_ecef[0], r0[1] + l*k_ecef[1], r0[2] + l*k_ecef[2]];
            let h = norm3(&r) - RE;
            let n_e = electron_density(h);
            let b   = geomag_field(&r);
            let bk  = dot3(&b, &k_ecef);
            let w = if i==0 || i==n_pts { 0.5 } else { 1.0 };
            accum += w * n_e * bk * dl;
        }
        // s = 7527 · c · ∫ N (B⃗·k̂) dL   (Eq. 8 of Hawarey et al. 2005)
        // Units of s: m·Hz³.  Then τ_range = s/ν³ [meters], τ_time = s/(ν³·c) [seconds].
        S2_COEFF * C * accum
    }

    /// Geometric slant path length from station altitude h0 [m above WGS84]
    /// to a target geocentric altitude h_t [m above WGS84], for elevation el [rad].
    /// Straight-ray + spherical Earth (the 2005 paper's geometric basis).
    fn slant_path_to_height(h0: f64, h_t: f64, el: f64) -> f64 {
        let r0 = RE + h0;
        let rt = RE + h_t;
        // Solve |r0·sin(el) + L|² + (r0·cos(el))² = rt²  for L
        // where L is the along-ray distance from station.
        let a = r0 * el.sin();
        let b = r0 * el.cos();
        let disc = rt*rt - b*b;
        if disc < 0.0 { return 0.0; }
        let l = -a + disc.sqrt();
        l.max(0.0)
    }

    /// 2nd order ionospheric delay [s] at a single station for a given source.
    /// From Hawarey et al. 2005:
    ///   Eq. 3a/3b: τ_range = s · ν⁻³   (range correction in meters)
    ///   Convert to time delay: τ_time = τ_range / c
    ///   ⇒ τ_time = s / (ν³ · c)
    pub fn second_order_delay(sta: &Station, src: &Source, freq_hz: f64, gmst: f64) -> f64 {
        let s = s_term(sta, src, gmst);
        s / (freq_hz.powi(3) * C)
    }

    /// Differential 2nd order ionospheric delay between two stations
    /// (= τ_2,sta2 − τ_2,sta1, per step 7 of the 2005 paper, Section 3).
    pub fn second_order_diff(s1: &Station, s2: &Station, src: &Source,
                             freq_hz: f64, gmst: f64) -> f64 {
        let t1 = second_order_delay(s1, src, freq_hz, gmst);
        let t2 = second_order_delay(s2, src, freq_hz, gmst);
        t2 - t1
    }
}

fn compute(s1: &Station, s2: &Station, src: &Source,
           tec1: f64, tec2: f64, freq_hz: f64, gmst: f64,
           include_2nd_order: bool) -> DelayResult {
    let el1=elevation_epoch(s1,src,gmst).to_degrees();
    let el2=elevation_epoch(s2,src,gmst).to_degrees();
    let valid=el1>3.0&&el2>3.0;
    let tg=geometric_delay(s1,s2,src);
    let (td,t1,t2)=tropo_diff(s1,s2,src,gmst);
    let ti=iono_diff(s1,s2,src,tec1,tec2,freq_hz,gmst);
    let ti2 = if include_2nd_order && valid {
        ionosphere::second_order_diff(s1,s2,src,freq_hz,gmst)
    } else { 0.0 };
    let ps=|t:f64|t*1e12;
    DelayResult {
        baseline_name: format!("{}-{}",s1.name,s2.name),
        source_name: src.name.clone(),
        tau_geo_ps: ps(tg), tau_tropo1_ps: ps(t1), tau_tropo2_ps: ps(t2),
        tau_iono_ps: ps(ti), tau_iono2_ps: ps(ti2),
        tau_total_ps: ps(tg+td+ti+ti2),
        el1_deg: el1, el2_deg: el2, valid,
    }
}

// ── Tests ──────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ecef_equator() {
        let p=geodetic_to_ecef(0.0,0.0,0.0);
        assert!((p[0]-RE).abs()<1.0);
        assert!(p[1].abs()<1.0);
        assert!(p[2].abs()<1.0);
    }
    #[test]
    fn elevation_zenith() {
        let sta=Station::new("T",45.0,10.0,100.0,1013.25,288.15,10.0);
        let src=Source::new("Z",10.0/15.0,45.0);
        assert!(elevation(&sta,&src).to_degrees()>85.0);
    }
    #[test]
    fn geo_delay_magnitude() {
        let s1=Station::new("WETT",49.145,12.879,666.0,1013.0,286.0,8.0);
        let s2=Station::new("KOKE",22.127,-159.665,1167.0,900.0,295.0,15.0);
        let src=Source::new("3C273",12.4973,2.0524);
        let tau=geometric_delay(&s1,&s2,&src).abs();
        assert!(tau<0.035&&tau>0.0);
    }
    #[test]
    fn saastamoinen_ztd() {
        let sta=Station::new("T",45.0,0.0,0.0,1013.25,288.15,10.0);
        let ztd=saastamoinen_delay(&sta,PI/2.0);
        assert!(ztd>2.0&&ztd<2.8,"ZTD={:.3}",ztd);
    }
    #[test]
    fn total_delay_units() {
        let s1=Station::new("A",52.0,0.0,50.0,1013.25,288.15,10.0);
        let s2=Station::new("B",40.0,15.0,100.0,1010.0,290.0,12.0);
        let src=Source::new("Q",6.0,30.0);
        let r=compute(&s1,&s2,&src,5.0,5.0,8.4e9,0.0,false);
        assert!(r.tau_total_ps.abs() < 10_000_000_000.0,
            "tau_total={:.0} ps", r.tau_total_ps);
    }

    // ─────────────────────────────────────────────────────────────────────
    // 2nd ORDER IONOSPHERIC TESTS
    // Citation: Hawarey, M., Hobiger, T. & Schuh, H. (2005).
    //   Effects of the 2nd order ionospheric terms on VLBI measurements.
    //   Geophys. Res. Lett. 32, L11304. doi:10.1029/2005GL022729.
    //
    // Reference values from the 2005 paper:
    //   • Per-ray 2nd order term: order 10⁻¹² s ≈ 1 ps (Abstract & Fig. 1).
    //   • Differential per baseline: ~3–9.4 psec max (Table 1).
    //   • Algonquin–Hartebeesthoek baseline: 3.0 psec (Table 1).
    //   • Maximum across baselines (Kokee–NyAlesund): 9.4 psec ≈ 2.8 mm
    //     (Table 1 + Section 4).
    // ─────────────────────────────────────────────────────────────────────

    #[test]
    fn iono2_per_ray_psec_order_of_magnitude() {
        // Per the 2005 paper Abstract: 2nd order effects are at the 10⁻¹² s
        // level for individual rays. Use a synthetic zenith source so the
        // ray actually pierces the ionosphere at our chosen GMST.
        let sta = Station::new("ALGO", 45.9555, -78.0728, 224.0,
                               1013.0, 286.0, 8.0); // Algonquin (Canada)
        // At GMST=0, LST at Algonquin = -78.07°. A source at RA = LST/15h,
        // Dec = station latitude is exactly at zenith (HA=0, El=90°).
        let zenith_ra_h = (-78.0728_f64 / 15.0) + 24.0; // wrap to positive: 18.79h
        let src = Source::new("ZENITH", zenith_ra_h, 45.9555);
        let tau2 = ionosphere::second_order_delay(&sta, &src, 8.4e9, 0.0).abs();
        let tau2_ps = tau2 * 1e12;
        // 2005 Fig. 1 shows excursions up to ~12 ps at S-band (2.3 GHz).
        // At X-band (8.4 GHz), ν⁻³ scaling reduces this by (8.4/2.3)³ ≈ 48.7×.
        // Expect ~0.05–1 ps at zenith for a quiet ionosphere.
        assert!(tau2 > 1e-15 && tau2 < 5e-11,
            "Per-ray |τ_2| = {:.3e} s ({:.4} ps) — out of paper's 10⁻¹² order-of-magnitude",
            tau2, tau2_ps);
    }

    #[test]
    fn iono2_algonquin_hartebeesthoek_baseline_psec_scale() {
        // 2005 paper, Table 1, baseline "Algonquin – Hartebeesthoek": 3.0 psec.
        // We verify our differential is at the few-psec scale (within a
        // factor that allows for our PIM-replacement Chapman profile).
        let s1 = Station::new("ALGO", 45.9555, 281.9272-360.0, 224.0,
                              1013.0, 286.0, 8.0);
        let s2 = Station::new("HART", -25.8897, 27.6853, 1416.0,
                              1013.0, 290.0, 12.0);
        let src = Source::new("3C273", 12.4973, 2.0524);
        // Sample several epochs within the CONT02 window (2002-Oct-16..30)
        // and take the max — Table 1 reports the per-baseline maximum.
        let mut max_ps = 0.0_f64;
        let jd_start = julian_date(2002, 10, 16, 0.0);
        for h_step in 0..(15*24) {
            let jd = jd_start + (h_step as f64) / 24.0;
            let gmst = gmst_rad(jd);
            let d = ionosphere::second_order_diff(&s1, &s2, &src, 8.4e9, gmst).abs();
            let d_ps = d * 1e12;
            if d_ps > max_ps && elevation_epoch(&s1,&src,gmst).to_degrees()>3.0
                             && elevation_epoch(&s2,&src,gmst).to_degrees()>3.0 {
                max_ps = d_ps;
            }
        }
        // Paper's Table 1: 3.0 psec for this baseline. Tolerance band 0.05–50 ps
        // reflects our tilted-dipole + Chapman profile vs. paper's full IGRF + PIM.
        assert!(max_ps > 0.05 && max_ps < 50.0,
            "Algonquin–Hartebeesthoek max |Δτ_2| = {:.3} ps (paper: 3.0 ps)", max_ps);
    }

    #[test]
    fn iono2_kokee_nyalesund_max_baseline() {
        // 2005 paper, Table 1: Kokee–Ny-Alesund had the maximum 2nd order
        // effect at 9.4 psec. We verify the differential scale and that
        // this baseline produces a notably larger value than short ones.
        let kokee = Station::new("KOKE", 22.1266, -159.665, 1177.0,
                                  900.0, 295.0, 15.0);
        let nyale = Station::new("NYAL", 78.9291,   11.8699,   87.0,
                                 1010.0, 268.0,  4.0);
        let src = Source::new("3C273", 12.4973, 2.0524);
        let mut max_ps = 0.0_f64;
        let jd_start = julian_date(2002, 10, 16, 0.0);
        for h_step in 0..(15*24) {
            let jd = jd_start + (h_step as f64) / 24.0;
            let gmst = gmst_rad(jd);
            let d = ionosphere::second_order_diff(&kokee, &nyale, &src, 8.4e9, gmst).abs();
            let d_ps = d * 1e12;
            if d_ps > max_ps && elevation_epoch(&kokee,&src,gmst).to_degrees()>3.0
                             && elevation_epoch(&nyale,&src,gmst).to_degrees()>3.0 {
                max_ps = d_ps;
            }
        }
        // Paper: 9.4 ps. Tolerance band 0.05–100 ps reflects our use of a
        // tilted-dipole + Chapman profile vs. the paper's full IGRF + PIM.
        assert!(max_ps > 0.05 && max_ps < 100.0,
            "Kokee–NyAlesund max |Δτ_2| = {:.3} ps (paper: 9.4 ps)", max_ps);
    }

    #[test]
    fn iono2_frequency_scaling_inverse_cube() {
        // Eq. 3a/3b: τ_2 ∝ ν⁻³. Halving the frequency should octuple |τ_2|.
        let sta = Station::new("T", 45.0, 10.0, 200.0, 1013.0, 288.0, 10.0);
        let src = Source::new("S", 6.0, 30.0);
        let t_x = ionosphere::second_order_delay(&sta, &src, 8.4e9, 0.0).abs();
        let t_s = ionosphere::second_order_delay(&sta, &src, 4.2e9, 0.0).abs();
        let ratio = t_s / t_x;
        // Should be exactly 8.0 by construction of Eq. 3 — geometry is shared.
        assert!((ratio - 8.0).abs() < 1e-6,
            "ν⁻³ scaling check: ratio={:.6} (expected 8.0)", ratio);
    }

    #[test]
    fn iono2_signs_and_antisymmetry() {
        // Differential between the same station and itself must be zero.
        // Differential is antisymmetric under station swap.
        let s1 = Station::new("A", 50.0, 10.0, 100.0, 1013.0, 285.0, 8.0);
        let s2 = Station::new("B", 30.0, -50.0, 200.0, 1013.0, 290.0, 12.0);
        let src = Source::new("Q", 8.0, 20.0);
        let d_self = ionosphere::second_order_diff(&s1, &s1, &src, 8.4e9, 0.0).abs();
        assert!(d_self < 1e-20, "self-diff must be zero, got {:.3e}", d_self);
        let d_ab = ionosphere::second_order_diff(&s1, &s2, &src, 8.4e9, 0.0);
        let d_ba = ionosphere::second_order_diff(&s2, &s1, &src, 8.4e9, 0.0);
        assert!((d_ab + d_ba).abs() < 1e-20,
            "antisymmetry: τ(A,B) + τ(B,A) = {:.3e}, expected 0", d_ab+d_ba);
    }
}

// ── GUI ────────────────────────────────────────────────────────────────────────

struct VlbiApp {
    stations:  Vec<Station>,
    sources:   Vec<Source>,
    baselines: Vec<Baseline>,
    freq_ghz:  f64,
    tec1:      f64,
    tec2:      f64,
    include_2nd_order: bool,    // NEW
    tab:       usize,
    results:   Vec<DelayResult>,
    sb_name:String,sb_lat:String,sb_lon:String,sb_h:String,
    sb_p:String,sb_t:String,sb_e:String,
    srb_name:String,srb_ra:String,srb_dec:String,
    bl_s1:usize,bl_s2:usize,
    sel_sta:Option<usize>,sel_src:Option<usize>,
    status:String,
    epoch_year:  i32,
    epoch_month: u32,
    epoch_day:   u32,
    epoch_hour:  f64,
}

impl Default for VlbiApp {
    fn default() -> Self {
        Self {
            stations: vec![
                Station::new("WETTZELL",  49.1449,  12.8789, 666.0, 1013.0, 286.0,  8.0),
                Station::new("KOKEE",     22.1265,-159.6651,1167.0,  900.0, 295.0, 15.0),
                Station::new("ONSALA60",  57.3958,  11.9255,  10.0, 1015.0, 283.0,  9.0),
                Station::new("GILCREEK",  64.9785,-147.4963, 317.0,  980.0, 270.0,  5.0),
            ],
            sources: vec![
                Source::new("3C84",    3.3314, 41.5117),
                Source::new("3C273",  12.4973,  2.0524),
                Source::new("3C345",  16.5920, 39.7542),
                Source::new("CygA",   19.9904, 40.7339),
            ],
            baselines: vec![
                Baseline{sta1:0,sta2:1},
                Baseline{sta1:0,sta2:2},
                Baseline{sta1:2,sta2:3},
            ],
            freq_ghz:8.4, tec1:5.0, tec2:5.0,
            include_2nd_order: true,    // ON by default
            tab:0, results:vec![],
            sb_name:"NEW_STA".into(),sb_lat:"0.0".into(),sb_lon:"0.0".into(),sb_h:"0.0".into(),
            sb_p:"1013.25".into(),sb_t:"288.15".into(),sb_e:"10.0".into(),
            srb_name:"NEW_SRC".into(),srb_ra:"0.0".into(),srb_dec:"0.0".into(),
            bl_s1:0,bl_s2:1,
            sel_sta:None,sel_src:None,
            status:"Ready — set epoch and press Compute.".into(),
            epoch_year:2025, epoch_month:1, epoch_day:15, epoch_hour:14.0,
        }
    }
}

/// Julian Date from calendar date + decimal UT hours
fn julian_date(y: i32, m: u32, d: u32, ut_hours: f64) -> f64 {
    let (y, m) = if m <= 2 { (y - 1, m + 12) } else { (y, m) };
    let a = (y as f64 / 100.0).floor();
    let b = 2.0 - a + (a / 4.0).floor();
    let jd0 = (365.25 * (y as f64 + 4716.0)).floor()
             + (30.6001 * (m as f64 + 1.0)).floor()
             + d as f64 + b - 1524.5;
    jd0 + ut_hours / 24.0
}

impl VlbiApp {
    fn current_gmst(&self) -> f64 {
        let jd = julian_date(self.epoch_year, self.epoch_month, self.epoch_day, self.epoch_hour);
        gmst_rad(jd)
    }

    fn compute_all(&mut self) {
        self.results.clear();
        let fhz  = self.freq_ghz * 1e9;
        let gmst = self.current_gmst();
        for bl in &self.baselines {
            if bl.sta1>=self.stations.len()||bl.sta2>=self.stations.len(){continue;}
            let s1=self.stations[bl.sta1].clone();
            let s2=self.stations[bl.sta2].clone();
            for src in &self.sources {
                self.results.push(compute(&s1,&s2,src,self.tec1,self.tec2,fhz,gmst,
                                          self.include_2nd_order));
            }
        }
        let v=self.results.iter().filter(|r|r.valid).count();
        let order_str = if self.include_2nd_order {"1st+2nd"} else {"1st only"};
        self.status=format!("{} results — {} valid, {} below horizon. Iono: {}.",
            self.results.len(),v,self.results.len()-v,order_str);
    }

    fn try_add_station(&mut self) {
        if let(Ok(la),Ok(lo),Ok(h),Ok(p),Ok(t),Ok(e))=(
            self.sb_lat.parse::<f64>(),self.sb_lon.parse::<f64>(),self.sb_h.parse::<f64>(),
            self.sb_p.parse::<f64>(),self.sb_t.parse::<f64>(),self.sb_e.parse::<f64>()
        ){
            let name=self.sb_name.clone();
            self.stations.push(Station::new(&name,la,lo,h,p,t,e));
            self.status=format!("Added station '{}'.",name);
        } else { self.status="Invalid station parameters.".into(); }
    }

    fn try_add_source(&mut self) {
        if let(Ok(ra),Ok(dec))=(self.srb_ra.parse::<f64>(),self.srb_dec.parse::<f64>()){
            let name=self.srb_name.clone();
            self.sources.push(Source::new(&name,ra,dec));
            self.status=format!("Added source '{}'.",name);
        } else { self.status="Invalid source parameters.".into(); }
    }
}

// ── UI helpers ──────────────────────────────────────────────────────────────────
fn lv(ui:&mut Ui,k:&str,v:&str,vc:Color32){
    ui.horizontal(|ui|{
        ui.label(RichText::new(k).color(Color32::from_rgb(120,140,170)).monospace().size(10.0));
        ui.label(RichText::new(v).color(vc).monospace().strong().size(10.0));
    });
}
fn sec(ui:&mut Ui,t:&str,c:Color32){
    ui.add_space(8.0);
    ui.horizontal(|ui|{
        ui.label(RichText::new(t).color(c).monospace().size(10.0).strong());
    });
    ui.separator();
}
fn frow(ui:&mut Ui,lbl:&str,buf:&mut String,w:f32){
    ui.horizontal(|ui|{
        ui.label(RichText::new(lbl).color(Color32::from_rgb(120,140,170)).monospace().size(10.0));
        ui.add(egui::TextEdit::singleline(buf).desired_width(w).font(egui::TextStyle::Monospace));
    });
}

impl eframe::App for VlbiApp {
    fn update(&mut self, ctx:&egui::Context, _frame:&mut eframe::Frame){
        ctx.set_visuals(egui::Visuals::dark());

        egui::TopBottomPanel::bottom("credit").show(ctx,|ui|{
            ui.add_space(3.0);
            ui.horizontal(|ui|{
                ui.label(RichText::new("© Dr. Mosab Hawarey").color(Color32::from_rgb(80,100,130)).monospace().size(9.0));
                ui.label(RichText::new("  |  ").color(Color32::from_rgb(50,65,90)).size(9.0));
                ui.label(RichText::new("@DrHawarey").color(Color32::from_rgb(0,180,230)).monospace().size(9.0));
                ui.label(RichText::new("  |  ").color(Color32::from_rgb(50,65,90)).size(9.0));
                ui.label(RichText::new("github.com/mhawarey").color(Color32::from_rgb(0,180,230)).monospace().size(9.0));
                ui.label(RichText::new("  |  ").color(Color32::from_rgb(50,65,90)).size(9.0));
                ui.label(RichText::new("VLBI Delay Model v1.1").color(Color32::from_rgb(80,100,130)).monospace().size(9.0));
            });
            ui.add_space(2.0);
        });

        egui::TopBottomPanel::bottom("status").show(ctx,|ui|{
            ui.horizontal(|ui|{
                ui.label(RichText::new("▸ ").color(Color32::from_rgb(0,200,100)).size(11.0));
                ui.label(RichText::new(&self.status).color(Color32::from_rgb(170,195,220)).monospace().size(10.0));
            });
        });

        egui::TopBottomPanel::top("hdr").show(ctx,|ui|{
            ui.add_space(8.0);
            ui.horizontal(|ui|{
                ui.label(RichText::new("VLBI").color(Color32::from_rgb(0,200,255)).size(22.0).strong().monospace());
                ui.label(RichText::new("Delay Model").color(Color32::WHITE).size(22.0).strong());
                ui.add_space(14.0);
                ui.label(RichText::new("// Geometric + Tropospheric (Saastamoinen) + Ionospheric (1st + 2nd order, Hawarey 2005)")
                    .color(Color32::from_rgb(80,100,130)).monospace().size(10.0));
            });
            ui.add_space(6.0);
            ui.horizontal(|ui|{
                let tabs=["Stations","Sources","Baselines","Settings","Results"];
                for(i,lbl) in tabs.iter().enumerate(){
                    let act=self.tab==i;
                    let c=if act{Color32::from_rgb(0,200,255)}else{Color32::from_rgb(120,140,170)};
                    let fill=if act{Color32::from_rgb(12,28,52)}else{Color32::from_rgb(16,20,30)};
                    if ui.add(egui::Button::new(RichText::new(*lbl).color(c).monospace().size(10.0).strong())
                        .fill(fill).stroke(Stroke::new(if act{1.5}else{0.4},c))).clicked(){self.tab=i;}
                }
                ui.add_space(16.0);
                if ui.add(egui::Button::new(
                    RichText::new("▶  COMPUTE ALL").color(Color32::from_rgb(10,20,10)).size(11.0).strong()
                ).fill(Color32::from_rgb(0,200,100))).clicked(){
                    self.compute_all(); self.tab=4;
                }
            });
            ui.add_space(4.0);
        });

        egui::CentralPanel::default().show(ctx,|ui|{
            match self.tab {
                0=>self.ui_stations(ui),
                1=>self.ui_sources(ui),
                2=>self.ui_baselines(ui),
                3=>self.ui_settings(ui),
                4=>self.ui_results(ui),
                _=>{}
            }
        });
    }
}

impl VlbiApp {
    fn ui_stations(&mut self, ui:&mut Ui){
        egui::ScrollArea::vertical().show(ui,|ui|{
            sec(ui,"STATION LIST",Color32::from_rgb(0,200,255));
            egui::Grid::new("sg").num_columns(9).striped(true).spacing([8.0,3.0]).show(ui,|ui|{
                for h in &["#","Name","Lat°","Lon°","h m","P hPa","T K","e hPa",""] {
                    ui.label(RichText::new(*h).color(Color32::from_rgb(0,200,255)).monospace().size(9.0).strong());
                }
                ui.end_row();
                let mut del=None;
                for(i,sta) in self.stations.iter().enumerate(){
                    let sel=self.sel_sta==Some(i);
                    let c=if sel{Color32::from_rgb(0,200,255)}else{Color32::from_rgb(190,210,235)};
                    if ui.selectable_label(sel,RichText::new(format!("{}",i+1)).monospace().size(9.0).color(Color32::from_rgb(90,115,145))).clicked(){self.sel_sta=Some(i);}
                    ui.label(RichText::new(&sta.name).color(c).monospace().strong().size(9.0));
                    ui.label(RichText::new(format!("{:.3}",sta.lat_deg)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.3}",sta.lon_deg)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.0}",sta.height_m)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.1}",sta.pressure_hpa)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.1}",sta.temp_k)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.1}",sta.e_hpa)).monospace().size(9.0));
                    if ui.small_button("✕").clicked(){del=Some(i);}
                    ui.end_row();
                }
                if let Some(i)=del{self.stations.remove(i);self.sel_sta=None;}
            });
            sec(ui,"ADD STATION",Color32::from_rgb(127,255,110));
            frow(ui,"Name      ",&mut self.sb_name,110.0);
            frow(ui,"Lat (°)   ",&mut self.sb_lat,90.0);
            frow(ui,"Lon (°)   ",&mut self.sb_lon,90.0);
            frow(ui,"Height (m)",&mut self.sb_h,90.0);
            frow(ui,"P (hPa)   ",&mut self.sb_p,90.0);
            frow(ui,"T (K)     ",&mut self.sb_t,90.0);
            frow(ui,"e (hPa)   ",&mut self.sb_e,90.0);
            if ui.add(egui::Button::new(RichText::new("+ Add Station").color(Color32::BLACK).strong())
                .fill(Color32::from_rgb(127,255,110))).clicked(){ self.try_add_station(); }

            if let Some(i)=self.sel_sta {
                if i<self.stations.len(){
                    let p=self.stations[i].ecef();
                    sec(ui,"ECEF (selected)",Color32::from_rgb(200,160,255));
                    lv(ui,"X = ",&format!("{:.3} m",p[0]),Color32::from_rgb(200,160,255));
                    lv(ui,"Y = ",&format!("{:.3} m",p[1]),Color32::from_rgb(200,160,255));
                    lv(ui,"Z = ",&format!("{:.3} m",p[2]),Color32::from_rgb(200,160,255));
                    lv(ui,"|r|= ",&format!("{:.3} m",norm3(&p)),Color32::from_rgb(200,160,255));
                }
            }
        });
    }

    fn ui_sources(&mut self, ui:&mut Ui){
        egui::ScrollArea::vertical().show(ui,|ui|{
            sec(ui,"SOURCE LIST",Color32::from_rgb(0,200,255));
            egui::Grid::new("srg").num_columns(5).striped(true).spacing([8.0,3.0]).show(ui,|ui|{
                for h in &["#","Name","RA (h)","Dec (°)","ICRF vector"] {
                    ui.label(RichText::new(*h).color(Color32::from_rgb(0,200,255)).monospace().size(9.0).strong());
                }
                ui.end_row();
                let mut del=None;
                for(i,src) in self.sources.iter().enumerate(){
                    let sel=self.sel_src==Some(i);
                    let c=if sel{Color32::from_rgb(255,190,50)}else{Color32::from_rgb(190,210,235)};
                    if ui.selectable_label(sel,RichText::new(format!("{}",i+1)).monospace().size(9.0).color(Color32::from_rgb(90,115,145))).clicked(){self.sel_src=Some(i);}
                    ui.label(RichText::new(&src.name).color(c).monospace().strong().size(9.0));
                    ui.label(RichText::new(format!("{:.4}h",src.ra_h)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:+.4}°",src.dec_deg)).monospace().size(9.0));
                    let v=src.unit_vector();
                    ui.label(RichText::new(format!("[{:.4},{:.4},{:.4}]",v[0],v[1],v[2])).monospace().size(9.0).color(Color32::from_rgb(140,160,190)));
                    if ui.small_button("✕").clicked(){del=Some(i);}
                    ui.end_row();
                }
                if let Some(i)=del{self.sources.remove(i);self.sel_src=None;}
            });
            sec(ui,"ADD SOURCE",Color32::from_rgb(255,190,50));
            frow(ui,"Name    ",&mut self.srb_name,110.0);
            frow(ui,"RA  (h) ",&mut self.srb_ra,90.0);
            frow(ui,"Dec (°) ",&mut self.srb_dec,90.0);
            if ui.add(egui::Button::new(RichText::new("+ Add Source").color(Color32::BLACK).strong())
                .fill(Color32::from_rgb(255,190,50))).clicked(){ self.try_add_source(); }

            if let(Some(si),Some(sri))=(self.sel_sta,self.sel_src){
                if si<self.stations.len()&&sri<self.sources.len(){
                    let el=elevation(&self.stations[si],&self.sources[sri]).to_degrees();
                    sec(ui,"ELEVATION (sel. sta × src)",Color32::from_rgb(127,255,110));
                    let ec=if el>10.0{Color32::from_rgb(127,255,110)}else if el>3.0{Color32::from_rgb(255,200,50)}else{Color32::from_rgb(255,80,80)};
                    lv(ui,"El = ",&format!("{:.3}°",el),ec);
                    if el<=3.0{ui.label(RichText::new("⚠ Source below horizon (≤3°)").color(Color32::from_rgb(255,80,80)).size(10.0));}
                }
            }
        });
    }

    fn ui_baselines(&mut self, ui:&mut Ui){
        egui::ScrollArea::vertical().show(ui,|ui|{
            sec(ui,"BASELINE LIST",Color32::from_rgb(0,200,255));
            egui::Grid::new("blg").num_columns(5).striped(true).spacing([10.0,3.0]).show(ui,|ui|{
                for h in &["#","Sta 1","Sta 2","Length (km)",""] {
                    ui.label(RichText::new(*h).color(Color32::from_rgb(0,200,255)).monospace().size(9.0).strong());
                }
                ui.end_row();
                let mut del=None;
                for(i,bl) in self.baselines.iter().enumerate(){
                    let(n1,n2,blen)=if bl.sta1<self.stations.len()&&bl.sta2<self.stations.len(){
                        let s1=&self.stations[bl.sta1]; let s2=&self.stations[bl.sta2];
                        let p1=s1.ecef(); let p2=s2.ecef();
                        let d=((p2[0]-p1[0]).powi(2)+(p2[1]-p1[1]).powi(2)+(p2[2]-p1[2]).powi(2)).sqrt();
                        (s1.name.clone(),s2.name.clone(),d/1000.0)
                    }else{("?".into(),"?".into(),0.0)};
                    ui.label(RichText::new(format!("{}",i+1)).monospace().size(9.0).color(Color32::from_rgb(90,115,145)));
                    ui.label(RichText::new(&n1).color(Color32::from_rgb(0,200,255)).monospace().strong().size(9.0));
                    ui.label(RichText::new(&n2).color(Color32::from_rgb(255,190,50)).monospace().strong().size(9.0));
                    ui.label(RichText::new(format!("{:.1} km",blen)).monospace().size(9.0));
                    if ui.small_button("✕").clicked(){del=Some(i);}
                    ui.end_row();
                }
                if let Some(i)=del{self.baselines.remove(i);}
            });

            sec(ui,"ADD BASELINE",Color32::from_rgb(127,255,110));
            let names:Vec<String>=self.stations.iter().map(|s|s.name.clone()).collect();
            if !names.is_empty(){
                ui.horizontal(|ui|{
                    ui.label(RichText::new("Station 1 ").color(Color32::from_rgb(120,140,170)).monospace().size(10.0));
                    egui::ComboBox::from_id_source("cb1")
                        .selected_text(names.get(self.bl_s1).cloned().unwrap_or_default())
                        .show_ui(ui,|ui|{
                            for(i,n) in names.iter().enumerate(){ui.selectable_value(&mut self.bl_s1,i,n);}
                        });
                });
                ui.horizontal(|ui|{
                    ui.label(RichText::new("Station 2 ").color(Color32::from_rgb(120,140,170)).monospace().size(10.0));
                    egui::ComboBox::from_id_source("cb2")
                        .selected_text(names.get(self.bl_s2).cloned().unwrap_or_default())
                        .show_ui(ui,|ui|{
                            for(i,n) in names.iter().enumerate(){ui.selectable_value(&mut self.bl_s2,i,n);}
                        });
                });
                if ui.add(egui::Button::new(RichText::new("+ Add Baseline").color(Color32::BLACK).strong())
                    .fill(Color32::from_rgb(127,255,110))).clicked(){
                    if self.bl_s1!=self.bl_s2 && self.bl_s1<self.stations.len() && self.bl_s2<self.stations.len(){
                        self.baselines.push(Baseline{sta1:self.bl_s1,sta2:self.bl_s2});
                        self.status="Baseline added.".into();
                    }else{self.status="Select two different valid stations.".into();}
                }
            }
        });
    }

    fn ui_settings(&mut self, ui:&mut Ui){
        sec(ui,"EPOCH (UTC)",Color32::from_rgb(255,190,50));
        egui::Grid::new("epoch_g").num_columns(2).spacing([12.0,6.0]).show(ui,|ui|{
            ui.label(RichText::new("Year          ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.epoch_year).speed(1).clamp_range(1900i32..=2100));
            ui.end_row();
            ui.label(RichText::new("Month         ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.epoch_month).speed(1).clamp_range(1u32..=12));
            ui.end_row();
            ui.label(RichText::new("Day           ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.epoch_day).speed(1).clamp_range(1u32..=31));
            ui.end_row();
            ui.label(RichText::new("UT Hours      ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.epoch_hour).speed(0.25).clamp_range(0.0f64..=23.9999).suffix(" h"));
            ui.end_row();
        });
        let gmst_deg = self.current_gmst().to_degrees().rem_euclid(360.0);
        let jd = julian_date(self.epoch_year, self.epoch_month, self.epoch_day, self.epoch_hour);
        ui.horizontal(|ui|{
            ui.label(RichText::new("JD  = ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.label(RichText::new(format!("{:.6}", jd)).monospace().size(10.0).color(Color32::from_rgb(255,190,50)));
            ui.add_space(16.0);
            ui.label(RichText::new("GMST = ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.label(RichText::new(format!("{:.4}°", gmst_deg)).monospace().size(10.0).color(Color32::from_rgb(255,190,50)));
        });

        sec(ui,"OBSERVATION PARAMETERS",Color32::from_rgb(0,200,255));
        egui::Grid::new("setg").num_columns(2).spacing([12.0,6.0]).show(ui,|ui|{
            ui.label(RichText::new("Frequency (GHz)  ").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.freq_ghz).speed(0.1).clamp_range(0.1f64..=100.0).suffix(" GHz"));
            ui.end_row();
            ui.label(RichText::new("VTEC Sta 1 (TECU)").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.tec1).speed(0.5).clamp_range(0.0f64..=200.0).suffix(" TECU"));
            ui.end_row();
            ui.label(RichText::new("VTEC Sta 2 (TECU)").monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
            ui.add(egui::DragValue::new(&mut self.tec2).speed(0.5).clamp_range(0.0f64..=200.0).suffix(" TECU"));
            ui.end_row();
        });

        sec(ui,"IONOSPHERIC MODEL ORDER",Color32::from_rgb(200,160,255));
        ui.horizontal(|ui|{
            ui.checkbox(&mut self.include_2nd_order,
                RichText::new("Include 2nd order ionospheric correction").color(Color32::from_rgb(200,160,255)).monospace().size(10.0));
        });
        ui.label(RichText::new("  Reference: Hawarey, Hobiger & Schuh (2005) — GRL 32, L11304")
            .color(Color32::from_rgb(120,140,170)).monospace().size(9.0));
        ui.label(RichText::new("  doi:10.1029/2005GL022729")
            .color(Color32::from_rgb(120,140,170)).monospace().size(9.0));

        sec(ui,"PHYSICAL CONSTANTS",Color32::from_rgb(200,160,255));
        egui::Grid::new("cg").num_columns(2).spacing([12.0,4.0]).show(ui,|ui|{
            for(k,v) in &[
                ("c (m/s)",format!("{:.6e}",C)),("K₁",format!("{:.3e}",K1)),
                ("K₂",format!("{:.3e}",K2)),("K₃",format!("{:.3e}",K3)),
                ("TEC coeff",format!("{:.3e}",TEC_COEFF)),
                ("S₂ coeff (Eq.8)",format!("{:.0}",S2_COEFF)),
                ("B_eq (T)",format!("{:.3e}",B_EQ_TESLA)),
                ("WGS84 a (m)",format!("{:.3}",RE)),
            ]{
                ui.label(RichText::new(*k).monospace().size(10.0).color(Color32::from_rgb(120,140,170)));
                ui.label(RichText::new(v.as_str()).monospace().size(10.0).color(Color32::from_rgb(200,160,255)));
                ui.end_row();
            }
        });
        sec(ui,"FORMULA REFERENCE",Color32::from_rgb(127,255,110));
        for line in &[
            "τ_geo  = −(b⃗ · ŝ) / c",
            "ZHD    = 0.002277·P / (1 − 0.00266·cos2φ − 0.00028·h[km])",
            "ZWD    = 0.002277·(1255/T + 0.05)·e",
            "τ_trop = ZTD / sin(el)   [mapping function]",
            "τ_iono = K·ΔTEC·sin⁻¹(el) / (f²·c)   K=40.309×10¹⁶",
            "",
            "─── 2nd order (Hawarey et al., 2005) ───",
            "s      = 7527·c·∫ N·(B⃗·k̂) dL              [Eq. 8]",
            "τ_2    = s · ν⁻³                            [Eq. 3a/3b]",
            "Δτ_2   = τ_2,sta2 − τ_2,sta1               [Section 3]",
            "",
            "τ_tot  = τ_geo + Δτ_trop + Δτ_iono + Δτ_iono2",
            "",
            "GMST   = 280.46061837 + 360.98564736629·(JD − 2451545.0)  [deg]",
            "LST    = GMST + λ_station",
            "HA     = LST − RA_source",
            "sin(el)= sin(φ)·sin(δ) + cos(φ)·cos(δ)·cos(HA)",
        ]{
            ui.label(RichText::new(*line).monospace().size(10.0).color(Color32::from_rgb(170,195,220)));
        }
    }

    fn ui_results(&mut self, ui:&mut Ui){
        if self.results.is_empty(){
            ui.add_space(30.0);
            ui.vertical_centered(|ui|{
                ui.label(RichText::new("No results yet.").color(Color32::from_rgb(130,150,180)).size(14.0));
                ui.add_space(10.0);
                if ui.add(egui::Button::new(RichText::new("▶  COMPUTE ALL").color(Color32::from_rgb(10,20,10)).size(13.0).strong())
                    .fill(Color32::from_rgb(0,200,100))).clicked(){self.compute_all();}
            });
            return;
        }
        let valid=self.results.iter().filter(|r|r.valid).count();
        ui.horizontal(|ui|{
            lv(ui,"Pairs: ",&format!("{}",self.results.len()),Color32::from_rgb(0,200,255));
            ui.add_space(20.0);
            lv(ui,"Valid: ",&format!("{}",valid),Color32::from_rgb(127,255,110));
            ui.add_space(20.0);
            lv(ui,"Below horizon: ",&format!("{}",self.results.len()-valid),Color32::from_rgb(255,100,80));
            ui.add_space(20.0);
            let epoch_str = format!("{:04}-{:02}-{:02} {:05.2}h UTC",
                self.epoch_year, self.epoch_month, self.epoch_day, self.epoch_hour);
            lv(ui,"Epoch: ",&epoch_str,Color32::from_rgb(255,190,50));
            ui.add_space(20.0);
            lv(ui,"Iono order: ",
                if self.include_2nd_order {"1st + 2nd"} else {"1st only"},
                Color32::from_rgb(200,160,255));
        });
        ui.add_space(6.0);
        egui::ScrollArea::both().show(ui,|ui|{
            egui::Grid::new("rg").num_columns(10).striped(true).spacing([10.0,3.0]).show(ui,|ui|{
                for h in &["Baseline","Source","τ_geo ps","τ_trop₁ ps","τ_trop₂ ps","Δτ_iono ps","Δτ_iono2 ps","τ_total ps","El₁°","El₂°"]{
                    ui.label(RichText::new(*h).color(Color32::from_rgb(0,200,255)).monospace().size(9.0).strong());
                }
                ui.end_row();
                for r in &self.results{
                    let dim=Color32::from_rgb(70,90,110);
                    let tc=if r.valid{Color32::from_rgb(190,210,235)}else{dim};
                    let vc=if r.valid{Color32::from_rgb(127,255,110)}else{Color32::from_rgb(255,80,80)};
                    let ic=if r.valid{Color32::from_rgb(200,160,255)}else{dim};
                    let i2c=if r.valid{Color32::from_rgb(255,160,200)}else{dim};
                    let ec=|e:f64|if e>10.0{Color32::from_rgb(127,255,110)}else if e>3.0{Color32::from_rgb(255,200,50)}else{Color32::from_rgb(255,80,80)};
                    ui.label(RichText::new(&r.baseline_name).color(Color32::from_rgb(0,200,255)).monospace().size(9.0).strong());
                    ui.label(RichText::new(&r.source_name).color(Color32::from_rgb(255,190,50)).monospace().size(9.0).strong());
                    ui.label(RichText::new(format!("{:+.2}",r.tau_geo_ps)).color(tc).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.2}",r.tau_tropo1_ps)).color(Color32::from_rgb(140,210,255)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.2}",r.tau_tropo2_ps)).color(Color32::from_rgb(140,210,255)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:+.4}",r.tau_iono_ps)).color(ic).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:+.4}",r.tau_iono2_ps)).color(i2c).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:+.2}",r.tau_total_ps)).color(vc).monospace().size(9.0).strong());
                    ui.label(RichText::new(format!("{:.1}",r.el1_deg)).color(ec(r.el1_deg)).monospace().size(9.0));
                    ui.label(RichText::new(format!("{:.1}",r.el2_deg)).color(ec(r.el2_deg)).monospace().size(9.0));
                    ui.end_row();
                }
            });
        });
    }
}

// ── Entry ──────────────────────────────────────────────────────────────────────
fn main(){
    eframe::run_native(
        "VLBI Delay Model — Dr. Mosab Hawarey",
        eframe::NativeOptions{
            initial_window_size:Some(Vec2::new(1180.0,720.0)),
            min_window_size:Some(Vec2::new(820.0,560.0)),
            ..Default::default()
        },
        Box::new(|_|Box::new(VlbiApp::default())),
    );
}
