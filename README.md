# VLBI Delay Model

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A native desktop GUI application for computing VLBI interferometric delays — geometric, tropospheric, and first-order ionospheric — across multiple baselines and radio sources simultaneously.**

**Built in Rust with [egui](https://github.com/emilk/egui) / [eframe](https://github.com/emilk/egui/tree/master/crates/eframe).**

## Screenshots

| Stations Tab | Sources Tab |
|---|---|
| ![Stations](vlbi_preview1.png) | ![Sources](vlbi_preview2.png) |

| Settings Tab (Epoch) | Results Tab |
|---|---|
| ![Settings](vlbi_preview3.png) | ![Results](vlbi_preview4.png) |

![Results Detail](vlbi_preview5.png)

## Delay Components

| Component | Model | Formula |
|---|---|---|
| **Geometric** | Exact VLBI formula | τ_geo = −(b⃗ · ŝ) / c |
| **Tropospheric** | Saastamoinen (1972) | ZTD / sin(el), with dry + wet |
| **Ionospheric** | First-order thin-shell | K · ΔTEC · sin⁻¹(el) / (f² · c) |
| **Total** | Sum | τ_total = τ_geo + Δτ_trop + Δτ_iono |

## Features

- **Batch computation** — all (baseline × source) pairs in one click
- **5 tabs**: Stations, Sources, Baselines, Settings, Results
- **Epoch-aware elevation** — GMST / LST / Hour Angle computed per station
- **Add/remove** stations and sources at runtime
- **ECEF preview** per station (WGS84)
- **Color-coded results** — green = valid, yellow = low elevation, red = below horizon
- **Settings** — drag-adjust epoch (year/month/day/UT), frequency, and TEC

## Tropospheric Model (Saastamoinen)

```
ZHD = 0.002277 · P / (1 − 0.00266·cos2φ − 0.00028·h[km])
ZWD = 0.002277 · (1255/T + 0.05) · e
ZTD = ZHD + ZWD
τ_trop = ZTD / sin(el)       [simple mapping function]
```

Required inputs per station: surface pressure P (hPa), temperature T (K), water vapour partial pressure e (hPa).

## Ionospheric Model (1st order)

```
τ_iono = K · ΔTEC / (f² · c)
ΔTEC = TEC₂/sin(el₂) − TEC₁/sin(el₁)   [thin-shell slant TEC]
K = 40.309 × 10¹⁶  [m·Hz²/TECU]
```

## Epoch & Elevation

```
GMST   = 280.46061837 + 360.98564736629 · (JD − 2451545.0)
LST    = GMST + λ_station
HA     = LST − RA_source
sin(el)= sin(φ)·sin(δ) + cos(φ)·cos(δ)·cos(HA)
```

## Default Stations

| Station | Lat (°) | Lon (°) | h (m) |
|---|---|---|---|
| WETTZELL | 49.145 | 12.879 | 666 |
| KOKEE | 22.127 | −159.665 | 1167 |
| ONSALA60 | 57.396 | 11.926 | 10 |
| GILCREEK | 64.978 | −147.496 | 317 |

## Build & Run

### Prerequisites
```
Rust ≥ 1.75 (stable)
```

On Linux, also install:
```bash
sudo apt install libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev \
                 libxkbcommon-dev libssl-dev pkg-config
```

### Build
```bash
cargo build --release
./target/release/vlbi_delay
```

On Windows, double-click `run.bat`.

### Tests
```bash
cargo test
# 5 passed; 0 failed
```

## References

- Sovers, O.J., Fanselow, J.L. & Jacobs, C.S. (1998). Astrometry and geodesy with radio interferometry. *Rev. Mod. Phys.* 70(4):1393.
- Saastamoinen, J. (1972). Atmospheric correction for the troposphere and stratosphere. *AGU Geophys. Monogr.* 15:247–251.
- Spilker, J.J. (1994). Tropospheric effects. In: *Global Positioning System: Theory and Applications*, Vol. 1.
- Bassiri, S. & Hajj, G.A. (1993). Higher-order ionospheric effects on the GPS observables. *Manuscripta Geodaetica* 18:280–289.

## Remark on 2nd-order ionospheric terms

- 2nd-order ionospheric terms are not taken into account in this software. Researchers who are interested in 2nd-order ionospheric terms are advised to consult the author's 2005 paper:
Hawarey, M., Hobiger, T., & Schuh, H. (2005). Effects of the 2nd order ionospheric terms on VLBI measurements. Geophysical Research Letters, 32(11), L11304. https://doi.org/10.1029/2005GL022729

## Author

**Dr. Mosab Hawarey**
>
PhD, Geodetic & Photogrammetric Engineering (ITU) | MSc, Geomatics (Purdue) | MBA (Wales) | BSc, MSc (METU)

- GitHub: https://github.com/mhawarey
- Personal: https://hawarey.org/mosab
- ORCID: https://orcid.org/0000-0001-7846-951X

## License

MIT License