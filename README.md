# Splines

Simple implementation of common 1D interpolating splines.

Included methods are linear, Hermite, Cardinal, finite difference, Catmull-Rom and Akima.

## Usage

```odin
import "splines"

centers := []f64{0.0, 0.24, 0.34, 0.64, 0.75}
values := []f64{0.1, 0.29, 0.38, 0.54, 0.85}

model := splines.build_linear(centers, values)
sample := splines.eval_linear(model, 0.45)
```
