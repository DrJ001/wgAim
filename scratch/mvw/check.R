
load("

## check mvwgaim against new code

base.qtl <- asreml(
    fixed     = yield ~ Trial + Trial:Rep,
    random    = ~ fa(Trial, 2):id,
    residual  = ~ dsum(~units | Trial),
    data      = pheno,
    trace     = FALSE,
    na.action = asreml::na.method(x = "include"),
    maxit     = 80)

diag.qtl <- asreml(
    fixed     = yield ~ Trial + Trial:Rep,
    random    = ~ diag(Trial):id,
    residual  = ~ dsum(~units | Trial),
    data      = pheno,
    trace     = FALSE,
    na.action = asreml::na.method(x = "include"),
    maxit     = 80)

class(genObj)[1] <- "interval"

mvw <- mvwgaim(diag.qtl, base.qtl, genObj, merge.by = "id", Trait = "Trial", gen.type = "interval", n.fa = 2,
breakout = 3)

args(mvwgaim.asreml)
