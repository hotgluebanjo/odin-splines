// Simple implementation of common 1D interpolating splines.
package splines

import "core:math"

Float :: f64

Linear :: struct {
    centers: []Float,
    values: []Float,
    extrapolate: bool,
}

build_linear :: proc(centers, values: []Float, extrapolate: bool = true) -> Linear {
    assert(len(centers) == len(values))
    assert(len(centers) >= 2)

    return Linear{centers, values, extrapolate}
}

eval_linear :: proc(s: ^Linear, x: Float) -> Float {
    oob_res, oob := handle_oob(s.centers, s.values, x, s.extrapolate)
    if oob {
        return oob_res
    }

    i := find_interval(s.centers, x)

    return norm_lerp(x, s.centers[i], s.values[i], s.centers[i+1], s.values[i+1])
}

Hermite :: struct {
    centers: []Float,
    values: []Float,
    tangents: []Float,
    extrapolate: bool,
}

// A Cubic Hermite Spline built with the auto-tangents `method`.
//
// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
build_hermite :: proc(centers, values: []Float, method: Hermite_Method, extrapolate: bool = true) -> Hermite {
    assert(len(centers) == len(values))
    assert(len(centers) >= 4)
    n := len(centers)

    tangents := make([]Float, n)

    // First and last points.
    tangents[0] = (values[1] - values[0]) / (centers[1] - centers[0])
    tangents[n-1] = (values[n-1] - values[n-2]) / (centers[n-1] - centers[n-2])

    switch method {
    case .Cardinal:
        for i in 1..<n-1 {
            tangents[i] = (values[i+1] - values[i-1]) / (centers[i+1] - centers[i-1])
        }
    case .Mean_Velocity:
        for i in 1..<n-1 {
            v_1 := (values[i] - values[i-1]) / (centers[i] - centers[i-1])
            v0 := (values[i+1] - values[i]) / (centers[i+1] - centers[i])
            tangents[i] = (v_1 + v0) / 2.0
        }
    case .Catmull_Rom:
        for i in 1..<n-1 {
            delta_1 := centers[i] - centers[i-1]
            delta0 := centers[i+1] - centers[i]

            v_1 := (values[i] - values[i-1]) / delta_1
            v0 := (values[i+1] - values[i]) / delta0

            tangents[i] = (delta0 * v_1 + delta_1 * v0) / (delta0 + delta_1)
        }
    case .Pchip:
        for i in 1..<n-1 {
            delta_1 := centers[i] - centers[i-1]
            delta0 := centers[i+1] - centers[i]

            v_1 := (values[i] - values[i-1]) / delta_1
            v0 := (values[i+1] - values[i]) / delta0

            wl := 2.0 * delta0 + delta_1
            wr := delta0 + 2.0 * delta_1

            tangents[i] = (wl + wr) / (wl / v_1 + wr / v0)
        }
    }

    return Hermite{centers, values, tangents, extrapolate}
}

eval_hermite :: proc(s: ^Hermite, x: Float) -> Float {
    oob_res, oob := handle_oob(s.centers, s.values, x, s.extrapolate)
    if oob {
        return oob_res
    }

    i := find_interval(s.centers, x)

    delta := s.centers[i+1] - s.centers[i]
    t := (x - s.centers[i]) / delta

    h00 := 2.0 * math.pow(t, 3.0) - 3.0 * math.pow(t, 2.0) + 1.0
    h10 := math.pow(t, 3.0) - 2.0 * math.pow(t, 2.0) + t
    h01 := -2.0 * math.pow(t, 3.0) + 3.0 * math.pow(t, 2.0)
    h11 := math.pow(t, 3.0) - math.pow(t, 2.0)

    return (
        h00 * s.values[i] +
        h10 * s.tangents[i] * delta +
        h01 * s.values[i+1] +
        h11 * s.tangents[i+1] * delta
    )
}

Hermite_Method :: enum {
    // Previous--next secant-line tangents. No tension parameter as it causes ripples
    // at any value other than `0.0` in 1D.
    //
    // https://www.youtube.com/watch?v=UCtmRJs726U
    // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Cardinal_spline
    Cardinal,

    // https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Finite_difference
    Mean_Velocity,

    // Correctly derived non-uniform Catmull--Rom tangents.
    //
    // https://splines.readthedocs.io/en/latest/euclidean/catmull-rom-properties.html
    Catmull_Rom,

    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html
    Pchip,
}

Akima :: struct {
    centers: []Float,
    values: []Float,
    coeff: [][4]Float,
    extrapolate: bool,
}

// https://en.wikipedia.org/wiki/Akima_spline
build_akima :: proc(centers, values: []Float, extrapolate: bool = true) -> Akima {
    assert(len(centers) == len(values))
    assert(len(centers) >= 5)
    n := len(centers)

    tangents := make([]Float, n-1)
    defer delete(tangents)

    for i in 0..<n-1 {
        tangents[i] = (values[i+1] - values[i]) / (centers[i+1] - centers[i])
    }

    // N slopes: one per point.
    slopes := make([]Float, n)
    defer delete(slopes)

    // First and last points.
    slopes[0] = tangents[0]
    slopes[n-1] = tangents[n-2]

    // Second and last intervals.
    slopes[1] = (tangents[0] + tangents[1]) / 2.0
    slopes[n-2] = (tangents[n-3] + tangents[n-2]) / 2.0

    for i in 2..<n-2 {
        mn := (
            abs(tangents[i+1] - tangents[i]) * tangents[i-1] +
            abs(tangents[i-1] - tangents[i-2]) * tangents[i]
        )
        md := abs(tangents[i+1] - tangents[i]) + abs(tangents[i-1] - tangents[i-2])

        if md == 0.0 {
            slopes[i] = (tangents[i-1] + tangents[i]) / 2.0
        } else {
            slopes[i] = mn / md
        }
    }

    coeff := make([][4]Float, n-1)

    for i in 0..<n-1 {
        coeff[i][0] = values[i]
        coeff[i][1] = slopes[i]
        coeff[i][2] = (3.0 * tangents[i] - 2.0 * slopes[i] - slopes[i+1]) / (centers[i+1] - centers[i])
        coeff[i][3] = (slopes[i] + slopes[i+1] - 2.0 * tangents[i]) / math.pow(centers[i+1] - centers[i], 2.0)
    }

    return Akima{
        centers,
        values,
        coeff,
        extrapolate,
    }
}

eval_akima :: proc(s: ^Akima, x: Float) -> Float {
    oob_res, oob := handle_oob(s.centers, s.values, x, s.extrapolate)
    if oob {
        return oob_res
    }

    i := find_interval(s.centers, x)

    a := s.coeff[i][0]
    b := s.coeff[i][1]
    c := s.coeff[i][2]
    d := s.coeff[i][3]
    dist := x - s.centers[i]

    return a + b * dist + c * math.pow(dist, 2.0) + d * math.pow(dist, 3.0)
}

// Boundary_Condition :: enum {
//     Natural,
//     Parabolic,
//     Secant,
//     Inner,
// }

// boundary_tangent :: proc(x0, x1, y0, y1, m: Float, c: Boundary_Condition) -> (res: Float) {
//     switch c {
//     case .Natural:
//         res = 3.0 * (y1 - y0) / (2.0 * (x1 - x0)) - m / 2.0
//     case .Parabolic:
//         res = 2.0 * (y1 - y0) / (x1 - x0) - m
//     case .Secant:
//         res = (y1 - y0) / (x1 - x0)
//     case .Inner:
//         res = m
//     }
//     return
// }

// Interval lower index: i for i <= x < i+1.
// Clamped to [0, n-2] inclusively.
find_interval :: proc(points: []Float, x: Float) -> int {
    lower, upper := 0, len(points) - 1
    for lower != upper - 1 {
        center := lower + (upper - lower) / 2
        if x >= points[center] {
            lower = center
        } else {
            upper = center
        }
    }
    return lower
}

// Lerp between `y0` and `y0` using `x` normalized from [x0, x1] to [0, 1].
norm_lerp :: proc(x, x0, y0, x1, y1: Float) -> Float {
    x_norm := (x - x0) / (x1 - x0)
    return math.lerp(y0, y1, x_norm)
}

// Handle out-of-bounds `x`.
handle_oob :: proc(centers, values: []Float, x: Float, extrapolate: bool) -> (res: Float, oob: bool) {
    n := len(centers)
    switch {
    case x <= centers[0]:
        if extrapolate {
            res = norm_lerp(x, centers[0], values[0], centers[1], values[1])
        } else {
            res = values[0]
        }
        oob = true
    case x >= centers[n-1]:
        if extrapolate {
            res = norm_lerp(x, centers[n-2], values[n-2], centers[n-1], values[n-1])
        } else {
            res = values[n-1]
        }
        oob = true
    case:
        oob = false
    }
    return
}
