// Simple implementation of common 1D interpolating splines.
package splines

import "core:math"

Float :: f64

Linear :: struct {
    centers: []Float,
    values: []Float,
    end_tangents: [2]Float,
    extrapolate: bool,
}

build_linear :: proc(centers, values: []Float, extrapolate: bool = true) -> Linear {
    assert(len(centers) == len(values))
    assert(len(centers) >= 2)
    n := len(centers)

    end_tangents := [2]Float{
        end_tangent(centers[0], centers[1], values[0], values[1], 0.0, .Slope),
        end_tangent(centers[n-2], centers[n-1], values[n-2], values[n-1], 0.0, .Slope),
    }

    return Linear{centers, values, end_tangents, extrapolate}
}

eval_linear :: proc(s: ^Linear, x: Float) -> Float {
    n := len(s.centers)
    beyond_val, is_beyond := handle_beyond_range(
        s.centers[0],
        s.centers[n-1],
        s.values[0],
        s.values[n-1],
        s.end_tangents[0],
        s.end_tangents[1],
        x,
        s.extrapolate,
    )

    if is_beyond {
        return beyond_val
    }

    i := find_interval(s.centers, x)
    t := (x - s.centers[i]) / (s.centers[i+1] - s.centers[i])

    return math.lerp(s.values[i], s.values[i+1], t)
}

Hermite :: struct {
    centers: []Float,
    values: []Float,
    coeff: [][4]Float,
    end_tangents: [2]Float,
    extrapolate: bool,
}

// A Cubic Hermite Spline built with the auto-tangents `method`.
//
// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
build_hermite :: proc(
    centers, values: []Float,
    method: Hermite_Method,
    ends: End_Condition = .Natural,
    extrapolate: bool = true,
) -> Hermite {
    assert(len(centers) == len(values))
    assert(len(centers) >= 4)
    n := len(centers)

    tangents := make([]Float, n)
    defer delete(tangents)

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
    case .Akima:
        weights := make([]Float, n-1)
        defer delete(weights)

        for i in 0..<n-1 {
            weights[i] = (values[i+1] - values[i]) / (centers[i+1] - centers[i])
        }

        // Second and last intervals.
        tangents[1] = (weights[0] + weights[1]) / 2.0
        tangents[n-2] = (weights[n-3] + weights[n-2]) / 2.0

        for i in 2..<n-2 {
            mn := (
                abs(weights[i+1] - weights[i]) * weights[i-1] +
                abs(weights[i-1] - weights[i-2]) * weights[i]
            )
            md := abs(weights[i+1] - weights[i]) + abs(weights[i-1] - weights[i-2])

            if md == 0.0 {
                tangents[i] = (weights[i-1] + weights[i]) / 2.0
            } else {
                tangents[i] = mn / md
            }
        }
    }

    // First and last points.
    end_tangents := [2]Float{
        end_tangent(centers[0], centers[1], values[0], values[1], tangents[1], ends),
        end_tangent(centers[n-2], centers[n-1], values[n-2], values[n-1], tangents[n-2], ends),
    }

    tangents[0] = end_tangents[0]
    tangents[n-1] = end_tangents[1]

    coeff := make([][4]Float, n-1)

    for i in 0..<n-1 {
        delta := centers[i+1] - centers[i]
        coeff[i][0] = (2.0 * values[i] + tangents[i] * delta - 2.0 * values[i+1] + delta * tangents[i+1])
        coeff[i][1] = (-3.0 * values[i] + 3.0 * values[i+1] - 2.0 * delta * tangents[i] - delta * tangents[i+1])
        coeff[i][2] = delta * tangents[i]
        coeff[i][3] = values[i]
    }

    return Hermite{centers, values, coeff, end_tangents, extrapolate}
}

eval_hermite :: proc(s: ^Hermite, x: Float) -> Float {
    n := len(s.centers)
    beyond_val, is_beyond := handle_beyond_range(
        s.centers[0],
        s.centers[n-1],
        s.values[0],
        s.values[n-1],
        s.end_tangents[0],
        s.end_tangents[1],
        x,
        s.extrapolate,
    )

    if is_beyond {
        return beyond_val
    }

    i := find_interval(s.centers, x)
    t := (x - s.centers[i]) / (s.centers[i+1] - s.centers[i])

    return s.coeff[i][0] * math.pow(t, 3.0) + s.coeff[i][1] * math.pow(t, 2.0) + s.coeff[i][2] * t + s.coeff[i][3]
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

    // https://en.wikipedia.org/wiki/Akima_spline
    Akima,
}

End_Condition :: enum {
    Natural,
    Parabolic,
    Slope,
    Inner,
}

// An end tangent found from the last two points and/or the innermost one's derivative, `m`.
end_tangent :: proc(x0, x1, y0, y1, m: Float, c: End_Condition) -> (res: Float) {
    switch c {
    case .Natural:
        res = 3.0 * (y1 - y0) / (2.0 * (x1 - x0)) - m / 2.0
    case .Parabolic:
        res = 2.0 * (y1 - y0) / (x1 - x0) - m
    case .Slope:
        res = (y1 - y0) / (x1 - x0)
    case .Inner:
        res = m
    }
    return
}

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

// Handle out-of-bounds `x`. Optionally extrapolates using the end derivatives
// `m0` and `mn_1`.
handle_beyond_range :: proc(x0, xn_1, y0, yn_1, m0, mn_1, x: Float, extrapolate: bool) -> (res: Float, oob: bool) {
    switch {
    case x <= x0:
        if extrapolate {
            res = (x - x0) * m0 + y0
        } else {
            res = y0
        }
        oob = true
    case x >= xn_1:
        if extrapolate {
            res = (x - xn_1) * mn_1 + yn_1
        } else {
            res = yn_1
        }
        oob = true
    }
    return
}
