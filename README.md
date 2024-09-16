# Splines

Simple implementation of common 1D interpolating splines.

Included methods are linear and several Hermite derivatives.

## Usage

```odin
import "splines"

centers := []f64{0.0, 0.24, 0.34, 0.64, 0.75}
values := []f64{0.1, 0.29, 0.38, 0.54, 0.85}

model_linear := splines.build_linear(centers, values)
sample_linear := splines.eval_linear(&model_linear, 0.45)

model_hermite := splines.build_hermite(centers, values, method = .Akima)
sample_hermite := splines.eval_hermite(&model_hermite, 0.45)
```
