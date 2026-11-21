# Pan-tissue single-cell loss of chromosome Y with detection-efficiency correction: pre-registered (21/11/2026, before any data look)

Data: CELLxGENE human datasets with normal male and female donors, 10x, >= 12 donors, adult ages. Pilot: kidney,
PBMC (199 donors), skin, skeletal muscle, white adipose, pancreas, bone marrow, retina/RPE, heart. Only cells with
disease = normal and raw integer counts.

Y set (male-specific region, broadly expressed): RPS4Y1 ZFY USP9Y DDX3Y UTY KDM5D EIF1AY TXLNGY NLGN4Y PRKY TMSB4Y.
Pseudo-Y set (negative control, true loss 0): autosomal genes matched gene by gene to the Y set on mean CP10k in
male cells (25 expression bins), 20 independent draws.
Per cell: total UMI N, Y UMI y, pseudo-Y UMI p.

Estimator (detection-efficiency unfolding): for each dataset x cell type x donor with >= 50 cells,
y ~ f * delta0 + (1 - f) * NB(mean mu N, size k), fitted by EM (k chosen by likelihood over a grid).
f = fraction of cells without Y. Naive comparator: fraction of cells with y = 0.
CI by 200 bootstrap resamples of cells.

Calibration (must pass, else NO-GO for the estimator):
C1 pseudo-Y: median f over male groups < 0.02 and 95th percentile < 0.05.
C2 female cells: median f for the Y set > 0.95 (ambient RNA may keep it below 1).
Hypotheses:
H1 (replication, blood): classical monocytes have the highest f among PBMC types with >= 50 donors,
and f rises with age in monocytes (donor-level slope > 0, permutation P < 0.05).
H2 (pan-tissue): among cell types with >= 10 male donors spanning >= 30 years, the age slope of f is positive
in more than half (sign test P < 0.05).
H3 (descriptive, new): rank cell types outside blood, brain and kidney by age slope; report those with slope
>= the monocyte slope.
Orthogonal check later: organ ranking vs YchrFISH scores (JCI Insight 2025, 11 organs).
Caveats fixed in advance: ambient Y RNA from intact cells biases f down; doublets bias f down; nuclei carry fewer Y
transcripts (lower mu), handled by the detection model; age is per donor and confounded with dataset.

Amendment (post hoc, after the first dataset, skin): male groups whose Y set is entirely zero are excluded
(Y genes unmeasured in that sub-study after integration), and groups with expected Y UMI per intact cell
mu * mean UMI < 2 are excluded as not identifiable. The same filter is applied to the pseudo-Y control.
Amendment 2 (post hoc, after skin and muscle): at expected Y UMI 2-5 the pseudo-Y control shows false zero inflation;
the reported quantity is the excess f(Y) - median f(pseudo-Y) of the same group, and the main threshold is expected Y >= 5.
