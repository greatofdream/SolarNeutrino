# Orbit
The existed work about the position of the sun
+ [`SolTrack`](https://han.vandersluys.nl/pub/PDFs/SolTrack.pdf): 
+ [`sunpos`](https://www.psa.es/sdg/sunpos.htm)
+ [`SPA`](https://doi.org/10.1016/j.solener.2003.12.003): package `sunposition`, used in `src/Orbit.py`
## Lagranian
+ Sun mass $M$, Earth mass $m$, radius $r$, angle $\theta$
$\mathcal{L}=\frac{1}{2}m(\dot{r}^2+r^2\dot{\theta}^2)+\frac{GMm}{r}$
+ $\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial\mathcal{L}}{\partial\dot{r}}-\frac{\partial\mathcal{L}}{\partial r}=0$: 
+ $\frac{\mathrm{d}}{\mathrm{d}t}\frac{\partial\mathcal{L}}{\partial\dot{\theta}}-\frac{\partial\mathcal{L}}{\partial \theta}=0$: $mr^2\frac{\mathrm{d}\dot{\theta}}{\mathrm{d}t}=0$
  + The angle speed is a constant, conservation of angular momentum .
+ The equaiton of motion: $r=\frac{ep}{1-e\cos\theta}$, 

