package test

import "core:fmt"

import "splines"

p4k_r := []f64{
    0.12651154398918152,
    0.17209523916244507,
    0.2372085005044937,
    0.3418610394001007,
    0.45308011770248413,
    0.5466615557670593,
    0.6243613958358765,
    0.6985563635826111,
    0.8392552137374878,
    0.9555448889732361,
}

p4k_g := []f64{
    0.123659148812294,
    0.16424179077148438,
    0.22333107888698578,
    0.3234575390815735,
    0.4240618646144867,
    0.5155408382415771,
    0.5889230966567993,
    0.6531014442443848,
    0.7934792041778564,
    0.9095908403396606,
}

p4k_b := []f64{
    0.12279334664344788,
    0.16511650383472443,
    0.22781753540039062,
    0.33257898688316345,
    0.439339280128479,
    0.533390998840332,
    0.6094257831573486,
    0.6798356771469116,
    0.8196576833724976,
    0.9356806874275208,
}

alexa_r := []f64{
    0.12316921353340149,
    0.16436529159545898,
    0.21895116567611694,
    0.3058355152606964,
    0.3883604109287262,
    0.4513855576515198,
    0.5025792717933655,
    0.5487018823623657,
    0.6219090223312378,
    0.6957394480705261,
}

alexa_g := []f64{
    0.12164027243852615,
    0.1614387333393097,
    0.21378368139266968,
    0.3018149137496948,
    0.3802224397659302,
    0.4445049464702606,
    0.49428239464759827,
    0.5397692918777466,
    0.6113743782043457,
    0.6840298175811768,
}

alexa_b := []f64{
    0.11265846341848373,
    0.1438184231519699,
    0.19718866050243378,
    0.2855287790298462,
    0.36468011140823364,
    0.4288386404514313,
    0.4788527488708496,
    0.5228954553604126,
    0.5948082804679871,
    0.6671851277351379,
}

main :: proc() {
    model_r := splines.build_cardinal(p4k_r, alexa_r)
    model_g := splines.build_cardinal(p4k_g, alexa_g)
    model_b := splines.build_cardinal(p4k_b, alexa_b)

    size := 500
    fmt.printf("LUT_1D_SIZE %d\n", size)
    fmt.printf("LUT_1D_INPUT_RANGE 0.0 1.0\n")
    for i in 0..=size {
        grid := f64(i) / f64(size)
        r := splines.eval_cardinal(&model_r, grid)
        g := splines.eval_cardinal(&model_g, grid)
        b := splines.eval_cardinal(&model_b, grid)
        fmt.printf("%0.8f %0.8f %0.8f\n", r, g, b)
    }
}
