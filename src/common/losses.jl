const h_robust=1.0

const c_cauchy=2.385 
const c_andrews = 1.339
const c_huber=1.345
const c_welsch=2.985
const c_fair = 1.4
const c_bisquare = 4.685
const c_logistic = 1.205
const c_talwar = 2.795

cauchy(r)   = 1 / (1 + r^2)
cauchy_sq(r) = 1 / (1 + r^2)^2
andrews(r)  = (abs(r)<pi) * sin(r) / r
huber(r)    = 1 / max(1.0, abs(r))
huber(δ,r) = 1 / max(δ, abs(r)) 
welsch(r)   = exp(-(r^2))
fair(r)     = 1 / (1 + abs(r))
bisquare(r) = (abs(r) < 1) * (1 - r^2)^2
logistic(r) = tanh(r)/r
talwar(r)   = 1 * (abs(r) < 1)