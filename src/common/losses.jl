const h_robust=1.0
const c_cauchy=2.385 
const c_huber=1.345
const c_welsch=2.985
const c_fair = 1.4

huber(r) = 1 / max(1.0, abs(r)) 
cauchy(r) = 1 / (1 + r^2)
welsch(r) = exp(-(r^2))
fair(r) = 1 / (1 + abs(r))