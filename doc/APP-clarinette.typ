#import "@local/modale:1.0.0": *

#show: modale.with(
  title: "APP Clarinette",
  subtitle: "Mécanique non linéaire de l'anche en-dedans : seuils d'oscillation et bifurcations",
  date: datetime(year: 2026, month: 11, day: 9),
  series: "APP",
  keywords: ("acoustique musicale", "non-linéarités", "oscillations", "bifurcations"),
  abstract: [
    On étudie la dynamique de l'anche d'une clarinette en régime basses fréquences, où la fréquence de jeu est très inférieure à la fréquence propre du premier mode d'anche. On étudie ensuite les régimes permanents et leur stabilité, ainsi que les bifurcations possibles.
  ],
  cover: true,
  cover-mode: 3,
  toc: true,
)

#symbols(
  ($y$, [Position de l'anche], [m]),
  ($omega_r$, [Pulsation propre], [rad/s]),
  ($q_r$, [Invers du facteur de qualité], [1]),
  ($f_s$, [Force sur l'anche], [N]),
  ($m_s$, [Masse modale], [kg]),
  ($K_s$, [Raideur intrinsèque], [N/m]),
  ($p_m$, [Pression dans la bouche], [Pa]),
  ($p$, [Pression dans le bec], [Pa]),
  ($Delta p$, [Différence de pression], [Pa]),
  ($H$, [Hauteur du bec], [m]),
  ($p_M$, [Pression de plaquage], [Pa]),
  ($S_j$, [Section du jet d'air], [m^2]),
  ($w$, [Largeur effective du canal d'anche], [m]),
  ($u$, [Débit d'écoulement dans le bec], [m^3/s]),
  ($u_r$, [Débit créé par la vibration de l'anche], [m^3/s]),
  ($S_r$, [Surface vibrante de l'anche], [m^2]),
  ($V_"eq"$, [Volume équivalent], [m^3]),
  ($Delta cal(l)$, [Correction de longueur acoustique], [m]),
  ($Z_c$, [Impédance caractéristique du tuyau], [Pa.s/m^3])
)
#pagebreak()

= Dynamique de l'anche
En ne considérant que le premier mode d'anche, on a un oscillateur amorti forcé par la pression en provenance de la bouche :
$
  y'' + q_r omega_r y' + omega_r^2 y = f_s / m_s
$
Avec
- $y$ la position du bout de l'anche
- $omega_r$ sa pulsation propre
- $q_r = 1/Q_r$ l'inverse du facteur de qualité
- $f_s = -Delta p = - (p_m - p)$ avec $p_m$ la pression dans la bouche et $p$ la pression dans le bec
- $m_s$ la masse modale du premier mode d'anche
- On note $K_s$ la raideur intrinsèque de l'anche, $K_s = m_s omega_r^2$

#hypothesis[
  On se place en régime basses fréquences, où la fréquence de jeu est très inférieure à la fréquence propre du premier mode d'anche.
]

Dans ce cadre, on peut négliger la dynamique de l'anche et la considérer comme un ressort statique :
$
  K_s y = -Delta p
$
On introduit la pression de plaquage $p_M$ définie par
$
  p_M = K_s H
$
Avec $H$ la hauteur du bec, telle que $p_M$ est la pression pour laquelle l'anche touche le bec ($y = -H$). On en déduit

#result(label: "Déplacement de l'anche")[
$
  y = -H (Delta p) / p_M
$ <mv_anche>
]

Ceci étant vrai tant que $y > -H$ (l'anche n'est pas plaquée contre le bec), ce qui est équivalent à $|Delta p| < p_M$.

= Écoulement --- effet valve
Le jet d'air provient de la bouche (où il est à la pression $p_m$ et à vitesse nulle) et s'écoule dans le canal d'anche. L'équation de Bernoulli permet d'écire la vitesse du jet en fonction de la pression dans le bec $p$, en considérant une séparation du jet en sortie du canal :
$
  v_j = sqrt((2 Delta p) / rho)
$
Avec $rho$ la masse volumique de l'air. En sortie du canal, on considère que le jet devient turbulent et se mélange. Son énergie cinétique est dissipée dans l'entrée du bec, d'où une source de débit par effet valve. On a la section du jet d'air :
$
  S_j = (y + H) w
$
Avec $w$ la largeur effective du canal d'anche. Le débit de l'écoulement dans le bec est donc donné par
$
  u = S_j v_j = cases(u_A (1 + y/H) sqrt((|Delta p|) / p_M) text("sgn")(Delta p) quad &text("si") y > -H,
  0 &text("si") y <= -H text("(anche battante)"))\
  "Avec" u_A = w H sqrt(2 / rho p_M)
$

= Débit d'anche et correction de longueur
Le jet d'air n'est pas la seule source de débit dans le bec, il faut aussi prendre en compte le débit créé par la vibration de l'anche : $u_r = -S_r (dif y)/(dif t)$ avec $S_r$ la surface vibrante de l'anche. Le débit total entrant dans le bec est donc
$
  u_("in") = u + u_r = [h star p]
$
Avec $h$ la transformée inverse de l'admittance d'entrée du bec.

#hypothesis[
  On considère que les paramètres du musicien sont constants pendant le jeu, notamment la pression dans la bouche $p_m$.
]

Sous les hypothèses basses-fréquences et $p_m$ constant, on peut écrire
$
  (dif y) / (dif t) = H / p_M (dif p) / (dif t)
$
On peut donc réécrire le débit d'anche :
$
  u_r = - V_"eq" / (rho c^2) (dif p) / (dif t) quad "où" V_"eq" = (rho c^2) / p_M S_r H
$

Et donc le débit total entrant dans le bec est
$
  u = u_"in" + V_"eq" / (rho c^2) (dif p) / (dif t)
$

On peut interpréter le terme $V_"eq" / (rho c^2)$ comme une compliance acoustique ajoutée à l'admittance d'entrée du bec. Cela peut se reformuler en termes de correction de longueur acoustique $V_"eq" = S Delta cal(l)$ avec $S$ la section du tuyau.

#result(label: "Correction de longueur du canal")[
  On a donc la correction de longueur effective du canal d'anche
  $
    Delta cal(l) = V_"eq" / S = (rho c^2) / p_M (S_r) / S H
  $
]

On en déduit donc l'expression du débit d'anche en fonction de la pression dans le bec et de la pression de plaquage :
$
  u = [h star p](t) + (S Delta cal(l)) / (H) p_M / (rho c^2) (dif y) / (dif t) quad "et" quad U(omega) = Y(omega) P(omega)
$

= Modèle global élémentaire à deux équations
Du fait de l'approximation de la dynamique de l'anche à celle d'un ressort statique, on peut supprimer une inconnue du problème (@mv_anche) et on arrive au système d'équations suivant :
$
  cases(
    u = [h star p](t) + (S Delta cal(l)) / (H) p_M / (rho c^2) (dif y) / (dif t),
    y = -H (Delta p) / p_M,
    u = cases(u_A (1 - (Delta p)/p_M) sqrt((|Delta p|) / p_M) text("sgn")(Delta p) quad &text("si") y > -H quad ("ou" Delta p <= p_M),
    0 &text("si") y <= -H quad ("ou" Delta p > p_M)) 
  )
$

On décide d'adimensionner ces équations afin d'alléger les notations. On introduit alors, avec $Z_c = rho c \/ S$ l'impédance caractéristique du tuyau, les variables  suivantes :
$
  tilde(y) = y / H = - (Delta p) / p_M; quad tilde(p) = p / p_M; quad tilde(u) = u Z_c / p_M; quad quad x = tilde(y) + gamma\
  gamma = p_m / p_M; quad zeta = Z_c u_A / p_M = Z_c w H sqrt(2 / (rho p_M)) = (w H) / S sqrt(2 C_p / C_v p_0 / p_M)
$

Où on a, dans l'expression de $zeta$, remplacé la compressibilité $1 \/ rho c^2$ par sa valeur en fonction de la pression atmosphérique $p_0$ et du rapport des capacités calorifiques $C_p \/ C_v$. Dans la suite on ignorera le tilde pour alléger.

Ici, $x$ correspond aux déplacement de l'anche par rapport à sa position d'équilibre sous pression de bouche constante $p_m$, et non par rapport à l'équilibre au repos. Le décalage $gamma$ correspond à l'enfoncement statique dû à $p_m$: en $p = 0$, on a $tilde(y) = -gamma$ et donc $x = 0$. De même, le débit $tilde(u)$ est le débit d'anche à l'équilibre sous condition de pression de bouche constante. Après subsitution et simplification, on se ramène aux équations suivantes :

$
  cases(
    u = [h star p](t) + (Delta cal(l)) / c (dif x) / (dif t),
    x = p,
    u = F(p) = cases(
      zeta (1 - gamma + p)sqrt(gamma - p) quad &"si" gamma - p <= 1,
      0 &"si" gamma - p >= 1
    )
  )
$

$p$ et $u$ sont la pression et le débit à l'entrée du résonateur (donc dans le bec). On peut encore un peu simplifier en intégrant la correction de longueur dans l'admittance d'entrée du bec $Y'(omega) = Y(omega) + j omega V_"eq" \/ (rho c^2)$ (valide seulement proche d'un pic d'impédance, à basses fréquences, pour un bec quasi-cylindrique, avec $K_s$ constant). On a donc

#result(label: "Système d'équations adimensionné")[
$
  cases(
    u = [h star p](t) quad "ou" quad U(omega) = Y'(omega) P(omega),
    u = F(p) = cases(
      zeta (1 - gamma + p)sqrt(gamma - p) quad &"si" gamma - p <= 1,
      0 &"si" gamma - p >= 1
    )
  )
$ <deux_eq>
]

Le problème est donc ramené à la résolution de la première équation, linéaire et avec mémoire, ainsi que $u = F(p)$, non-linéaire et ponctuelle.

= Régime statique --- méthode de résolution
Dans le cas où les inconnues $p$ et $u$ ne dépendent pas du temps, on a une solution triviale au système d'équations précédent. Si on suppose nulle l'impédance d'entrée à fréquence nulle (lorsque le débit est constant), on a également $p = 0$ et
$
  F_0 = F(0) = zeta (1 - gamma)sqrt(gamma)
$
Il n'y a pas de son (et dans le cas où l'anche bat, si $gamma >= 1$, $F_0 = 0$). Cependant, il existe deux cas de figure où le système présente un transitoire, même avec $u$ constant :
- Si la pression dans le bec avant qu'on souffle n'est pas nulle
- Si $gamma$ ne passe pas instantanément de $0$ à une valeur stable (ce qui est le cas réel)

Fondamentalement, si on part d'un état initial hors équilibre, même avec une excitation constante, le système présente une évolution vers un état stable.

On décide de réécrire le problème à l'aide des ondes de pression aller et retour dans le bec :
$
  p = p^+ + p^-; quad u = p^+ - p^-; quad p^+ = 1/2 (p + u); quad p^- = 1/2 (p - u)
$

L'équation non-linéaire $u = F(p)$ devient :
$
  p^+ - p^- = F(p^+ + p^-)
$

Et on suppose qu'on peut réécrire cette équation sous la forme
$
  p^+ = G(p^-)
$
De sorte que $G$ est la fonction non linéaire qui donne l'onde aller (injectée dans le conduit du bec) en fonction de l'onde retour (changée de signe et réfléchie en sortie du bec). Le signe $-$ vient du fait que $p^- = -p^+(t-tau)$ sans pertes ($p^-$ est l'ancienne onde aller retardée).

En réalité, il suffit d'effectuer une rotation de $45$° (et de multiplier par un facteur $1\/sqrt(2)$) pour passer de la représentation $(p,u)$ à la représentation $(p^+,p^-)$ (et donc d'obtenir $G$ à partir de $F$).

Avec ce changement de variables, on peut écrire la première équation du système :
$
  p^-(t) = [r star p^+](t) quad "ou" quad P^-(omega) = R(omega) P^+(omega)
$

Avec $R(omega)$ le coefficient de réflexion et $r(t)$ la fonction de réflexion. L'intérêt principal est que la fonction de réflexion est beaucoup plus courte que la réponse impulsionnelle et peut donc être employée pour un calcul numérique.

On peut en déduire l'équation
#result[
$
  p(t) = u(t) + [r star (p + u)](t) = u(t) + integral_0^t r(tau) [p(t-tau) + u(t-tau)]dif tau
$ <p_u>
]

En y injectant la seconde équation du système (@deux_eq), on en déduit $p$ et $u$ en fonction des conditions initiales. On simplifie dans la suite le problème pour étudier les mécanismes de production du son : on néglige les pertes et on ne considère qu'un seul mode du résonateur.

= Approximation sans pertes pour un cyclindre
Si on néglige l'effet des pertes et du rayonnement, la fonction de réflexion est simplement
$
  r(t) = -delta(t - tau), quad tau = (2 cal(l)) / c
$
On en déduit $p^- (t) = -p^+(t - tau)$ et on peut écrire @p_u comme
$
  p(t) = u(t) - p(t - tau) - u(t - tau)
$ <recurrence>

En discrétisant le temps de sorte que $p_n = p(n tau)$, on obtient le système
$
  cases(
    p_n^- = -p_(n-1)^+,
    p_n^+ = G(-p_n^-) = G(p_(n-1)^+),
  )
$
Avec la condition initiale $p_0^- = 0$, puisqu'aucune onde ne peut revenir plus vite que $c$. Les conditions initiales dictent l'allure de $G$ et donc l'évolution du système vers un état stable. Plusieurs régimes sont possibles.

== Stabilité du régime statique

On considère le cas où l'attaque est instantanée, c'est-à-dire que $gamma$ passe de $0$ à sa valeur finale en un instant. Au début du transitoire, $p$ est très petit, donc on peut linéariser l'équation $u = F(p)$ discrétisée autour de $0$ :
$
  u_n = F(p_n) approx F(0) + F'(0) p_n = F_0 + A p_n quad "où" quad A = zeta (3 gamma - 1) / (2 sqrt(gamma))  "et" F_0 = zeta (1 - gamma)sqrt(gamma)
$

Avec la relation $p(t) = u(t) - p(t - tau) - u(t - tau)$, on obtient
$
  p_n = -p_(n-1) (1 + A) / (1 - A) = p_0 (-1)^n ((1 + A) / (1 - A))^n quad "et" quad p_0 = p_0^+ = F_0 / (1 - A)
$

On remarque que le terme élevé à la puissance $n$ tend vers $0$ si $A$ est négatif (on tend vers le régime statique stable) et grandit exponentiellement si $A$ est positif (le régime statique est instable). On en déduit que $A = 0$ est le seuil d'oscillation. D'après l'expression de $A$, cela correspond à $gamma = 1\/3$, c'est-à-dire à une pression de bouche égale à un tiers de la pression de plaquage.

#result(label: "Condition d'oscillation")[
  Sous les hypothèses précédentes, le seuil d'oscillation est atteint lorsque la pression de bouche est égale à un tiers de la pression de plaquage : $p_m = p_M \/ 3$.
]

En supposant que l'anche ne bat pas et en poursuivant le développement autour de $0$ de $u = F(p)$, on a
$
  u = F(p) approx F_0 + A p + B p^2 + C p^3 \
  F_0 = zeta (1 - gamma)sqrt(gamma); quad A = zeta (3 gamma - 1) / (2 sqrt(gamma)); quad B = -zeta (3 gamma + 1) / (8 gamma^(3\/2)); quad C = - zeta (gamma + 1) / (16 gamma^(5\/2))
$

En première approximation, on remarque que tous les coefficients sont proportionnels à l'ouverture de l'anche au repos $zeta$, que $A$ varie quasi-linéairement avec $gamma$ et que $B$ et $C$ sont quasi-constants (dans la zone de variation de $gamma$).

On s'intéresse maintenant à l'existence d'un régime permanent périodique. On peut supposer en premier lieu sa période égale à $tau = 2 cal(l) \/ c$. On avait déjà (sans pertes) que $p^+(t - tau) = -p^-(t) = p^+(t)$ avec hypothèse de périodicité. On en déduit que $p(t) = p^+(t) + p^-(t) = 0$ et il n'y aurait pas de solution périodique non nulle.

En supposant une période de longueur $2tau$ on obtient $p(t) = -p(t+tau)$ (un signal à deux états dont les deux valeurs sont de signes opposés) et par conséquence $u(t) = u(t+tau)$ (le débit est donc constant en régime permanent). Les points limites $plus.minus p_infinity$ sont tels que $F(p_infinity) = F(-p_infinity)$ i.e.
$
  p_infinity = cases(
    plus.minus [(3 gamma - 1) (1 - gamma)]^(1\/2) &" si" gamma < 1\/2,
    plus.minus gamma &" si" gamma > 1\/2
  )
$

Dans notre cas $C$ est négatif, donc cette solution existe pour $A$ positif (donc lorsque le régime statique est instable).

== Stabilité du régime permanent
On étudie l'impact d'une perturbation autour de $p_(infinity)$ :
$
  u_n = F(p_n) approx F(p_(infinity)) + (p_n - p_(infinity)) F'(p_(infinity)) \
  u_(n+1) = F(p_(n+1)) approx F(-p_(infinity)) + (p_(n+1) + p_(infinity)) F'(-p_(infinity))
$
On note $F'^+ = F'(p_(infinity))$ et $F'^- = F'(-p_(infinity))$.
En utilisant l'équation $p_n + p_(n+1) = u_n - u_(n+1)$, on a
$
  (1 - F'^-)(p_(n+1) + p_infinity) = -(p_n - p_infinity) (1 + F'^+)\
  (1 - F'^+)(p_(n+2) - p_infinity) = -(p_(n+1) + p_infinity) (1 + F'^-)
$
D'où
$
  (p_(n+2) - p_infinity) = (p_n - p_infinity) ((1 + F'^+) (1 + F'^-)) / ((1 - F'^+) (1 - F'^-))
$
Pour qu'il y ait stabilité, il faut que le module du facteur multiplicatif soit inférieur à $1$ (l'amplitude de la perturbation diminue à chaque cycle). Ceci est équivalent à
$
  (F'^+ + F'^-)(1 + F'^+ F'^-) < 0
$
En reprenant l'expression du DL à l'ordre 3, on a
$
  A [1 + 4A(B^2/C + A)] > 0
$

On remarque que pour $A$ petit, la condition de stabilité est envron $A > 0$ et est opposée à la condition de stabilité du régime statique ($A < 0$).

Quand $A$ grandit, le terme du crochet linéaire en $A$ est $4A B^2 \/ C$, dont le signe dépend de celui de $C$:
- Si $C < 0$ (notre cas), le terme linéaire décroît lorsque $A$ augmente. En particulier, le crochet peut devenir négatif pour une valeur positive de $A$. Le régime permanent n'est plus stable. L'important c'est qu'on passe continûment du régime statique ($A = 0$) au régime oscillant ($A > 0$) en augmentant $A$ aussi lentement qu'on le veut (c'est une bifurcation dite "directe"). On peut donc jouer pianissimo sur une clarinette.
- Si $C > 0$, le régime permanent n'existe que pour $A < 0$, i.e. quand le régime statique est stable. Si $A$ est négatif et petit, le terme du crochet est positif et le régime permanenent est instable. Cependant il existe des valeurs de $A$ négatives plus grandes pour lesquelles la condition de stabilité du régime permanent est vérifiée. On a donc une bifurcation "inverse".

== Doublement de période & chaos
On peut montrer que le seuil d'existence d'une solution périodique de période $4 tau$ est bien
$
  A [ 1 + 4A(B^2/C + A)] = 0
$
On a donc une perte de stabilité du régime carré à ce moment là et l'apparition d'un régime de double période. Pour le modèle initial sans développement limité, on peut montrer qu'il existe des régimes de fréquence divisée par $2$, $4$, $6$, $12$, etc., dont les seuils en $gamma$ dépendent de $zeta$. En prenant en compte les pertes, certains doublements de période disparaissent. Enfin, en imposant une annulation de l'onde de pression en $t = tau \/ 3$ et une valeur non-nulle en $t = 0$ et $t = 2 tau \/ 3$ (en pratique en ouvrant un trou au tiers de la longueur), on peut obtenir des régimes de fréquence divisée par $3$, $6$, $9$, etc.

= Approximation à un mode de résonateur
En prenant en compte les pertes (sans le rayonnement) dans le résonateur et en tronquant l'impédance d'entrée à un seul mode, on peut écrire
$
  j omega F_1 U(omega) = (omega_1^2 - omega^2 + j omega omega_1 / Q_1) P(omega)
$
Avec $F_1 = 2 c \/ cal(l)$. Dans le domaine temporel, on a donc
$
  (dif^2 p) / (dif t^2) + omega_1 / Q_1 (dif p) / (dif t) + omega_1^2 p = F_1 (dif u) / (dif t)
$

En reprenant le développement à l'ordre de $3$ de $u$ autour de $p = 0$, on obtient
$
  (dif^2 p)/(dif t^2) + F_1 (dif) / (dif t) [(Y_(m 1) - A)p(t) - B p^2 (t) - C p^3 (t)] + omega_1^2 p(t) = 0
$

Où $Y_(m 1) = omega_1 \/ (F_1 Q_1)$ est la valeur du minimum réel de l'admittance d'entrée du résonateur.

== Stabilité du régime statique & début du transitoire
Une solution triviale est obtenue pour $p(t) = 0$, auquel cas, pour étudier sa stabilité (les petites variations autour de $0$), on linéarise l'équation en ignorant les termes d'ordre supérieur ($B$ et $C$) :
$
  (dif^2 p)/(dif t^2) + F_1 (Y_(m 1) - A) (dif p)/(dif t) + omega_1^2 p = 0
$

Il s'agit de l'équation d'un oscillateur amorti dont l'amortissement peut être négatif si $A > Y_(m 1)$ auquel cas les oscillations grandissent exponentiellement (la saturation de l'amplitude est alors assurée par les termes d'ordre supérieur). Dans le cas opposé, l'oscillation est exponentiellement armortie et le système est ramené au régime statique.

Il s'agit de la même condition de stabilité du régime statique que dans le cas sans résonateur, mais avec un seuil d'oscillation à $A = Y_(m 1)$ au lieu de $A = 0$ (il faut souffler plus fort).

== Existence du régime périodique
On suppose l'existence d'une solution sinusoïdale (pas rigoureux à cause de la présence d'une non-linéarité, mais fonctionne très bien en pratique). On applique la méthode de l'équilibrage harmonique pour obtenir l'amplitude de la solution sinusoïdale et on obtient :
$
  p(t) = 2 sqrt(1/3 (Y_(m 1) - A) / C) sin(omega_1 t + phi)
$

la condition d'existence est la même que la condition d'instabilité du régime statique, à savoir $A > Y_(m 1)$ (puisque $C$ est négatif). 