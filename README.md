[![ci-python-unittest](https://github.com/joebarteam11/EGR0D/actions/workflows/main.yml/badge.svg)](https://github.com/joebarteam11/EGR0D/actions/workflows/main.yml)

# Description

Based on a given configuration (3 TPX tanks [fuel, oxidizer and EGR],
a specific equivalence ratio, a percentage [mass, molar, volumetric]
of EGR and a kinetic scheme), the tool calculates equilibrium, 0D
reactors with specified residence time, or 1D premixed flames using
CANTERA, depending on the user's requirements.

Ranges of pressures, temperatures, richness and %EGR can be quickly
addressed since this tool can be run using mpi4py. Output data are
stored in Pandas dataframes and saved for post-processing.

# Definitions

Here are some key definitions used in this tool.

## Equivalence ratio

The equivalence ratio is defined as
$$\phi=S_{x}\frac{X_{fuel}}{X_{O_{2}}}$$

For Methane,
$$S_{x}=\left.\frac{X_{O_{2}}}{X_{fuel}}\right|_{stoech}=2$$

(hardcoded value for now, enhancement would be at list to be able to
provide it with other parameters)

## EGR rate

For validation purposes (see references section below), the EGR rate has
been defined as

$$\%EGR=X_{EGR}^{fuel}=\frac{n_{EGR}}{n_{fuel}+n_{EGR}}$$

but one parameter ('egr_unit') allows the user to select between 'mole',
'mass' or 'vol' definition. Moreover, the "egr_def" parameter allows to
choose between these definitions : TODO

## Flame thickness

The flame thickness is definied as the thermal thickness :
$$\delta_{T}=\frac{T_{out}-T_{in}}{\max\left(\overrightarrow{\nabla}T\right)}$$

# How to use

# Future developments

-   Add another fuel tank for CH4/H2 blend control loop

# Theoretical elements

## Enthalpy balance (CO2 reinjection only):

Consider the reaction equation (at stoichiometry):
$$CH_{4}+2(O_{2}+3,76N_{2})+zCO_{2}\rightarrow\left(1+z\right)CO_{2}+7,52N_{2}+2H_{2}O$$

Enthalpy balance :
$$H_{prod}\left(T_{P}\right)=H_{in}\left(T_{in}\right)=H_{reac}\left(T_{R}\right)+H_{CO_{2},r\acute{e}inj}\left(T_{P}\right)\Longleftrightarrow\underset{H_{prod-CO_{2}}\left(T_{P}\right)}{\underbrace{H_{prod}\left(T_{P}\right)-H_{CO_{2},r\acute{e}inj}\left(T_{P}\right)}}=H_{reac}\left(T_{R}\right)$$

With:
$`
\begin{cases} \ H_{reac}\left(T\right)=1h_{CH_{4}}^{m}\left(T\right)+2h_{O_{2}}^{m}\left(T\right)+7,52h_{N_{2}}^{m}\left(T\right) \ \\ \ H_{CO_{2},r\acute{e}inj}\left(T\right)=zh_{CO_{2}}^{m}\left(T\right) \ \\ \ H_{prod,\,sans\,pr\acute{e}l\grave{e}v}\left(T\right)=1h_{CO_{2}}^{m}\left(T\right)+2h_{H_{2}O}^{m}\left(T\right)+7,52h_{N_{2}}^{m}\left(T\right) \ \end{cases}
`$


and $H_{prod-CO_{2}}\left(T\right)=\left(1-z\right)h_{CO_{2}}^{m}\left(T\right)+2h_{H_{2}O}^{m}\left(T\right)+7,52h_{N_{2}}^{m}\left(T\right)$\
\
where $h^{m}\left(T\right)$ are the molar enthalpies of the species
calculated from NASA polynomials.\
Finally, $T_{P}$ is determined by linear interpolation between two
temperature limits $T_{a}$ and $T_{b}$ such that : 
$` \begin{aligned}T_{p} & =T_{a}+\frac{H_{prod-CO_{2}}\left(T_{P}\right)-H_{prod-CO_{2}}\left(T_{a}\right)}{H_{prod-CO_{2}}\left(T_{b}\right)-H_{prod-CO_{2}}\left(T_{a}\right)}\left(T_{b}-T_{a}\right)\\\Longleftrightarrow & T_{p}=T_{a}+\frac{H_{reac}\left(T_{R}\right)-H_{prod-CO_{2}}\left(T_{a}\right)}{H_{prod-CO_{2}}\left(T_{b}\right)-H_{prod-CO_{2}}\left(T_{a}\right)}\left(T_{b}-T_{a}\right)\end{aligned}  `$

## Calculating T from LHV without EGR

The two regimes are distinguished:

Rich :
$$\phi CH_{4}+2(O_{2}+3,76N_{2})\rightarrow CO_{2}+2H_{2}O+7,52N_{2}+\left(\phi-1\right)CH_{4}$$

Lean :
$$\phi CH_{4}+2(O_{2}+3,76N_{2})\rightarrow CO_{2}+2\phi H_{2}O+7,52N_{2}+\left(1-\phi\right)O_{2}$$

The stoichiometric mass ratio $s_{Y}$:
$$s_{Y}=s_{X}\frac{W_{air}}{W_{fuel}}=17.16$$

Calculation of LHV:

$$LHV\left[kJ/kg\right]=\frac{-\Delta H_{reaction}^{0}}{W_{fuel}}=\frac{-\left(H_{f,prod}^{0}-H_{f,reac}^{0}\right)\left[kJ/mol\right]}{W_{CH_{4}}\left[kg/mol\right]}$$

Calculation of flue gas $\overline{C_{p,GB}}\left(\phi\right)$:

$$\overline{C_{p,GB}}\left(\phi\right)=\frac{\underset{k}{\sum}\nu_{k}C_{p,mol,k}}{\underset{k}{\sum}\nu_{k}W_{k}}$$

Calculation of T :
$$T_{P}=T_{R}+\frac{LHV}{\overline{C_{p,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{tot}}$$

with
$\dot{m}\_{tot}=\dot{m}\_{f}+\dot{m}\_{air}=\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)$

In rich : $\dot{m}\_{f,br\hat{u}l\acute{e}}=\frac{\dot{m}\_{f}}{\phi}$

In lean: $\dot{m}\_{f,br\hat{u}l\acute{e}}=\dot{m}\_{f}$

## Calculation of T with EGR at $T_{p}$

Only the expression of $\dot{m}\_{tot}$ is modified:
$\dot{m}\_{tot}=\dot{m}\_{f}+\dot{m}\_{air}+\dot{m}\_{GB}$ (meaning
$\dot{m}\_{GB,r\acute{e}inj}$)

and $\dot{m}\_{GB}=\left(\dot{m}\_{f}+\dot{m}\_{air}\right)\frac{p}{1-p}$
with p the mass percentage of burnt gas reinjected, i.e. $p=\frac{\dot{m}\_{GB}}{\dot{m}\_{tot}}$ (see the section [below](#reasoning-on-the-molar-quantities-of-gas-reinjected-for-a-fixed-egr-rate) for a different reasoning).

Thus, $\dot{m}\_{tot}$ expression is : 
$' \begin{aligned}
\dot{m}\_{tot} & =\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)+\dot{m}\_{GB}=\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)+\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{p}{1-p}\right)\\
 & =\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(1+\frac{p}{1-p}\right)=\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)
\end{aligned} '$

$T_{p}$ expression become :
$$T_{P}=T_{mix}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{tot}}=T_{mel}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}$$

Where $T_{mix}$ is defined from:

$$T_{mix}\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]=T_{R}\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+T_{p}p\overline{C_{p,GB}}\left(\phi\right)$$

$$T_{P}=T_{mix}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}=\frac{T_{R}\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+T_{p}p\overline{C_{p,GB}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}$$

$$\Longleftrightarrow T_{P}-\frac{T_{p}p\overline{C_{p,GB}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}=\frac{T_{R}\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}$$

$$\Longleftrightarrow T_{P}=\frac{\frac{T_{R}\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}+\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}}{1-\frac{p\overline{C_{p,GB}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}}$$

$$T_{P}=T_{R}+\frac{\frac{LHV}{\overline{C_{P,GB}}\left(\phi\right)}\frac{\dot{m}\_{f,br\hat{u}l\acute{e}}}{\dot{m}\_{f}\left(1+\frac{s_{Y}}{\phi}\right)\left(\frac{1}{1-p}\right)}}{1-\frac{p\overline{C_{p,GB}}\left(\phi\right)}{\left[\left(1-p\right)\overline{C_{p,GF}}\left(\phi\right)+p\overline{C_{p,GB}}\left(\phi\right)\right]}}$$

## Reasoning on the molar quantities of gas reinjected for a fixed EGR rate

Let p be the [molar] percentage of sample in the burnt gases
reinjected at the inlet.

Let the reaction equation be the initial state (without re-injection):

$$CH_{4}+2(O_{2}+3,76N_{2})+0CO_{2}+0N_{2}\rightarrow\underset{\overbrace{z_{0}}}{1}CO_{2}+7,52N_{2}+2H_{2}O$$

Reasoning by recurrence :

$$CH_{4}\,+\,2(O_{2}+3,76N_{2})\,+\underset{\overbrace{p\times\left(1+z_{0}\right)}}{z_{1}}CO_{2}\,+\,z_{1}N_{2}\rightarrow(1+z_{1})CO_{2}\,+\,(1+z_{1})\times7,52N_{2}\,+\,2H_{2}O$$

$$CH_{4}\,+\,2(O_{2}+3,76N_{2})\,+\,z_{n}CO_{2}\,+\,z_{n}N_{2}\rightarrow(1\,+z_{n})CO_{2}\,+\,(1+z_{n})\times7,52N_{2}\,+\,2H_{2}O$$

This leads to the following relationship:

$$z_{n+1}=p\left(1+z_{n}\right)=pz_{n}+p$$

Which is the expression of an arithmetic-geometric sequence. Since
$p\in\left[0,1\right]$ then there exists a fixed point $\alpha$ such
that $\underset{n\rightarrow\infty}{\lim}z_{n}=\alpha$ with
$f(\alpha)=\alpha$ and $f(x)\hookrightarrow px+p$

Finally
$$\underset{n\rightarrow\infty}{\lim}z_{n}=\alpha=\frac{p}{1-p}$$

# References
test
