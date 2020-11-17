Example data and code for the manuscript "Flexible fitting of PROTAC concentration-response curves with Gaussian Processes". The kernel implemented in the Stan file, is a changepoint kernel 

$f1(x) ∼ GP(0, k1)$,  $f2(x) ∼ GP(0, k2)$, $f(x) := (1 − w(x))f1(x) + w(x)f2(x) ∼ GP(0, k)$, where $k(x, x')) = (1 − w(x))k1(x, x')(1 − w(x')) + w(x)k2(x, x)w(x')$ and we use a steep sigmoidal as the transition function: $w(x) = \sigma(kx), k=10$.
