# Theory

Here The basic theory of UPG and AAPDG will be explained.

## 1.UPG

UPG basically use a linear gravity assumption, which makes it easy to compute final state.

### 1.1 State equations

The equations of the system with linear gravity could be written as:
$$
\dot{\boldsymbol{r}} = \boldsymbol{v}\\
\dot{\boldsymbol{v}} = A\boldsymbol{r} + \boldsymbol{a_T}\\
\dot{m} = -\frac{T}{v_{ex}}
$$
where $A=-\mu/R_0^3$is the approximated gravity coefficient, $\boldsymbol{a_T}$ is the command acceleration, $T$ is the thrust, and $v_{ex}$ is the exhaust velocity.
The magnitude of $\boldsymbol{r}$ and $\boldsymbol{v}$ are very large, which can lead to diffiuclties when solving equations. A good practice is scale the postion down by $R_0$, velocity by $\sqrt{R_0g_0}$, and time by $\sqrt{\frac{R_0}{g_0}}$, where $R_0$ is the radius of the celestial body, and $g_0$ is the surface gravity. If you take derivatives of the scaled variables against scaled time, the state equations will be:
$$
\dot{\boldsymbol{r}} = \boldsymbol{v}\\
\dot{\boldsymbol{v}} = -\boldsymbol{r} + \frac{T}{m(t)g_0}\boldsymbol{1}_T\\
\dot{m} = -\frac{T}{v_{ex}}
$$
where $\boldsymbol{1}_T$ is the unit thrust vector.

### 1.2 The landing problems

The powered descent goal can be written as such
$$
\boldsymbol{r}_f = \boldsymbol{r}_{target}\\
\boldsymbol{v}_f = \boldsymbol{v}_{target}\\
$$
and we want to use as little fuels as possible, that is, to minimize
$$
J_1 = \int_0^{t_f} \frac{T}{v_{ex}} \mathrm{d}t
$$
This is called the pinpoint landing problem. However, this target can be difficult to solve, and therefore we can use a less strict contrain, which only requires the craft to land somewhere near the target. The new contraints are:
$$
\|\boldsymbol{r}_f\| = \|\boldsymbol{r}_{target}\|\\
\boldsymbol{v}_f = \boldsymbol{v}_{target}\\
$$
and to minimize
$$
J_2 = \kappa (\boldsymbol{r}_f - \boldsymbol{r}_{target})^T(\boldsymbol{r}_f - \boldsymbol{r}_{target}) + \int_0^{t_f} \frac{T}{v_{ex}} \mathrm{d}t
$$
where $\kappa > 0$ is a coefficient representing howimportant the landing difference is. Because this kind of performance index is called Bolza problem in optimal control problem, this problem is call bolza landing problem. It is obvious that the large $\kappa$ is, the closer the landing will be, and the more fuel will likely to be used for corrections. It is proven that the bolza landing always use less fuel than the pinpoint one. If $\kappa$ is set to 0, this yields a soft landing problem where the landing precision is not considered at all.

### 1.3 The solution for final states

To solver the above mentioned landing problems, one needs to calculate the final state of the craft when time reaches $t_f$.

It is known that with the optimal control theory, we have Hamiltonian $H$:
$$
H=\boldsymbol{p}_r^T \boldsymbol{v}+
\boldsymbol{p}_v^T \left( -\boldsymbol{r} + \frac{T}{m g_0}\boldsymbol{1}_T \right)-
p_m\frac{T}{v_{ex}} \sqrt{\frac{R_0}{g_0}}-
\frac{T}{v_{ex}} \sqrt{\frac{R_0}{g_0}}
$$
where $\boldsymbol{p}_r \in R^3$, $\boldsymbol{p}_v \in R^3$ and $p_m \in R$ are the costate vectors associated with $\boldsymbol{r}$, $\boldsymbol{v}$ and $m$, respectively. Accroding to optimal control theory, the thrust vector's direction should maximize $H$, which leads to
$$\boldsymbol{1}_T = \frac{\boldsymbol{p}_v}{\|\boldsymbol{p}_v\|}$$
and the costate vectors $\boldsymbol{p}_r$ and $\boldsymbol{p}_v$ satisfies
$$
\begin{pmatrix}
\dot{\boldsymbol{p}_r}\\
\dot{\boldsymbol{p}_v}
\end{pmatrix}=-
\begin{pmatrix}
\frac{\partial H}{\partial \boldsymbol{r}}\\[0.4em]
\frac{\partial H}{\partial \boldsymbol{v}}\\
\end{pmatrix}=
\begin{pmatrix}
\boldsymbol{p}_v\\
-\boldsymbol{p}_r
\end{pmatrix}
$$
Noticing that the costate vectors have the same equations as the state$\begin{pmatrix}
\dot{\boldsymbol{r}}\\
\dot{\boldsymbol{v}}
\end{pmatrix}$, we can immediately write the state transition matrix
$$
\boldsymbol{\Phi}(t-t_0)=
\begin{pmatrix}
\cos(t-t_0) \boldsymbol{I}_{3\times3} & \sin(t-t_0) \boldsymbol{I}_{3\times3}\\
-\sin(t-t_0) \boldsymbol{I}_{3\times3} & \cos(t-t_0) \boldsymbol{I}_{3\times3}
\end{pmatrix}
$$
where $\boldsymbol{I}_{3\times3}$ is the 3×3 identity matrix. With the transition matrix, we can write the final state at any time:
$$
\boldsymbol{x}=\boldsymbol{\Phi}(t-t_0) \boldsymbol{x}_0+
\int_{t_0}^t \boldsymbol{\Phi}(t-t_0) \boldsymbol{B} \frac{T}{m(t)g_0}\boldsymbol{1}_T \mathrm{d}t
$$
and the costate vectors with the same form:
$$
\boldsymbol{\lambda}=\boldsymbol{\Phi}(t-t_0) \boldsymbol{\lambda}_0
$$
where 
$\boldsymbol{x}=\begin{pmatrix} \boldsymbol{r}\\\boldsymbol{v} \end {pmatrix}$ is the state vector, 
$\boldsymbol{\lambda}=\begin{pmatrix} \boldsymbol{p}_r\\\boldsymbol{p}_v \end {pmatrix}$
is the costate vectors, 
$\boldsymbol{B}=\begin{pmatrix} 0\\1 \end{pmatrix}$
is the input matrix. To solve the integral, define
$$
\boldsymbol{I}_c(t,t_0)=\int_{t_0}^t \boldsymbol{1}_T(\tau)\cos(\tau)\frac{T}{m(\tau)g_0}\mathrm{d}\tau\\
\boldsymbol{I}_s(t,t_0)=\int_{t_0}^t \boldsymbol{1}_T(\tau)\sin(\tau)\frac{T}{m(\tau)g_0}\mathrm{d}\tau\\
$$
and
$$\boldsymbol{\Gamma}(t)=
\begin{pmatrix}
\sin(t-t_0) \boldsymbol{I}_{3\times3} & -\cos(t-t_0) \boldsymbol{I}_{3\times3}\\
\cos(t-t_0) \boldsymbol{I}_{3\times3} & \sin(t-t_0) \boldsymbol{I}_{3\times3}
\end{pmatrix}
$$
we can simplify the expression to
$$
\boldsymbol{x}=\boldsymbol{\Phi}(t-t_0) \boldsymbol{x}_0+
\boldsymbol{\Gamma}(t)\boldsymbol{I}(t,t_0)
$$
where
$\boldsymbol{I}(t,t_0)=
\begin{pmatrix}
\boldsymbol{I}_c(t,t_0)\\
\boldsymbol{I}_s(t,t_0)
\end{pmatrix}
$ is called thrust integrals, and can be numerically evaluated. For example, the Milne’s rule gives
$$
\boldsymbol{I}_c(t,t_0)=\frac{t-t_0}{90}
\left[7\boldsymbol{i}_c(t_0)+
32\boldsymbol{i}_c(t_0 + \delta)+
12\boldsymbol{i}_c(t_0 + 2\delta)+
32\boldsymbol{i}_c(t_0 + 3\delta)+
7\boldsymbol{i}_c(t_0 + 4\delta)
\right]
$$
where $\boldsymbol{i}_c$ is the integrand of $\boldsymbol{I}_c$, and $\delta=\frac{t-t_0}{4}$ is the time interval. For $\boldsymbol{I}_s(t,t_0)$, similar results can be obtained. Hence we can analytically calculate the state at any given time. Noticing that we have different phases with different thrust, the states should be calculated separately.

### 1.4 Constrains for root finding problem

To solve the landing problems, we need to find the root of the system, with 7 variables $\boldsymbol{z}=
\begin{pmatrix}
\boldsymbol{p}_r\\
\boldsymbol{p}_v\\
t_f
\end{pmatrix}$. The target contraints have given 6 equations for pinpoint landing and 4 constraints for bolza landing, so we still need to find 1 more constraint for pinpoint and 3 more for bolza.

The seventh constraint comes from the transversality conditions, which is $H(t_f)=0$ for free-time problem. Also, as $m(t_f)$ is not constrained, we have $p_m(t_f)=0$. This gives
$$
H(t_f)=\boldsymbol{p}_r^T(t_f) \boldsymbol{v}(t_f)-
\boldsymbol{p}_v^T(t_f) \boldsymbol{r}(t_f) + \frac{T}{m g_0}\|\boldsymbol{p}_v(t_f)\| -
\frac{T}{v_{ex}} \sqrt{\frac{R_0}{g_0}} = 0
$$
However, in this expression, the scales of $\boldsymbol{p}_r$ and $\boldsymbol{p}_v$ are very large, so we can use a equivilent version instead:
$$
H(t_f)=\boldsymbol{p}_r^T(t_f) \boldsymbol{v}(t_f)-
\boldsymbol{p}_v^T(t_f) \boldsymbol{r}(t_f) + \frac{T}{m g_0}\|\boldsymbol{p}_v(t_f)\| - 1= 0
$$
This equation plus the final postion and velocity constraints constitute a root finding problem for pinpoint landing.

For bolza landing, the last 2 equations comes from the transversality conditions associated with the terminal constraints,
$$
\boldsymbol{\lambda}(t_f)=-\kappa\frac{\partial \phi}{\partial \boldsymbol{x}_f}+\left[\frac{\partial \boldsymbol{s(\boldsymbol{x}_f)}}{\partial \boldsymbol{x}_f}\right]^T\boldsymbol{\nu}
$$
where $\phi = [\boldsymbol{r}(t_f)-\boldsymbol{r}_{target}]^T[\boldsymbol{r}(t_f)-\boldsymbol{r}_{target}]$, $\boldsymbol{s(\boldsymbol{x}_f)}$ is the 4 constraints of final state in bolza landing problem, and $\boldsymbol{\nu} \in R^4$ is a contant multiplier. To eliminate $\boldsymbol{\nu}$, one can arrive at the so called reduced transversality conditions:
$$
\boldsymbol{y}_i^T(\boldsymbol{x}_f)
\left(
\boldsymbol{\lambda}(t_f)+\kappa\frac{\partial \phi}{\partial \boldsymbol{x}_f}
\right)=0, i=1,2
$$
where $\boldsymbol{y}_i(\boldsymbol{x}_f) \in R^6$ are linearly independent solutions for the following equations:
$$
\frac{\partial \boldsymbol{s(\boldsymbol{x}_f)}}{\partial \boldsymbol{x}_f}\boldsymbol{y}=0
$$
assume that the third component of $\boldsymbol{r}(t_f)$ is not 0, we can obtain 2 solutions:
$$
\boldsymbol{y}_1=\begin{pmatrix}
1\\0\\-r_1(t_f)/r_3(t_f)\\0\\0\\0
\end{pmatrix},
\boldsymbol{y}_2=\begin{pmatrix}
0\\1\\-r_2(t_f)/r_3(t_f)\\0\\0\\0
\end{pmatrix}
$$
and the corresponding reduced transversality conditions are:
$$
r_3(t_f)p_{r_1}(t_f)-r_1(t_f)p_{r3}(t_f)+\kappa(r_1(t_f)r_{target_3}-r_3(t_f)r_{target_1})=0\\
r_3(t_f)p_{r_2}(t_f)-r_2(t_f)p_{r3}(t_f)+\kappa(r_2(t_f)r_{target_3}-r_3(t_f)r_{target_2})=0
$$
when $r_3(t_f)=0$, similar results can be obtained.

### 1.5 Jacobians

To solve the root finding problem, one can use the well developed algorithms, such as powell's method, lm, etc. It is found that powell's method works well. Note that the final states are analytically calculated, one can write the jacobians with close form to save time for numerically evaluate jacobians.
In this part, the Jacobians of the target function are given.

#### **1.5.1 Summary of the equations**

For convenience, we will breifly summarize how the target function computed.\

Given the target state $\boldsymbol{r}^*$ and $\boldsymbol{v}^*$, we should first calculate the final state $\boldsymbol{x}$ corresponds to the unknowns $\boldsymbol{z} = (\boldsymbol{p}_v^T \; \boldsymbol{p}_r^T \; t_f)^T \in R^7$. As the thrust will be 3-arc, we will have to calculate the states one by one.

Firstly, we calculate the state at time $t_1$:
$$
\boldsymbol{x}_1 =
\boldsymbol{\Phi}(t_1 - t_0) \boldsymbol{x}_0 +
\boldsymbol{\Gamma}(t_1) \, \boldsymbol{I}(t_1,t_0)
$$

Then we do the same thing for $t_2$ and $t_f$:

$$
\boldsymbol{x}_2 =
\boldsymbol{\Phi}(t_2 - t_1) \boldsymbol{x}_1 +
\boldsymbol{\Gamma}(t_2) \, \boldsymbol{I}(t_2,t_1)\\

\boldsymbol{x}_f =
\boldsymbol{\Phi}(t_f - t_2) \boldsymbol{x}_2 +
\boldsymbol{\Gamma}(t_f) \, \boldsymbol{I}(t_f,t_2)
$$

After we get the final state, we can calculate the 7 constraints:

$$
s_1=\boldsymbol{r}_f^T \boldsymbol{r}_f -
\boldsymbol{r}^{*T} \boldsymbol{r}^*\\[0.4em]

s_2=\boldsymbol{v}_{f1} - \boldsymbol{v}^*_1\\[0.4em]

s_3=\boldsymbol{v}_{f2} - \boldsymbol{v}^*_2\\[0.4em]

s_4=\boldsymbol{v}_{f3} - \boldsymbol{v}^*_3\\[0.4em]

s_5=\boldsymbol{p}_{rf}^T  \boldsymbol{v}_f -

\boldsymbol{p}_{vf}^T \boldsymbol{r}_f +

\frac{T_{max}}{m_f \, g_0} \| \boldsymbol{p}_{vf} \| -\boldsymbol{1}\\[0.4em]

s_6= \boldsymbol{y}_1^T \left[ \boldsymbol{p}_{rf} +
2\kappa (\boldsymbol{r}
_f-\boldsymbol{r}^*) \right]\\[0.4em]

s_7=\boldsymbol{y}_2^T \left[ \boldsymbol{p}_{rf} +
2\kappa (\boldsymbol{r}_f-\boldsymbol{r}^*) \right]\\
$$
here the * represent the target positon and velocity.

#### **1.5.2 Jacobians of thrust integral**

Now that we have obtained the target function $s=(s_1, s_2, s_3, s_4, s_5, s_6, s_7)^T$, we have to calculate the Jacobian
$$
\boldsymbol{J}=\frac{\partial \boldsymbol{s}}{\partial \boldsymbol{z}}
$$
If we derivate $x$ against $\lambda$, we can write the results one by one:
$$
\frac{\partial \boldsymbol{x}_1}{\partial \boldsymbol{\lambda}}
=\boldsymbol{\Gamma}(t_1)
\frac{\partial \boldsymbol{I}_1}{\partial \boldsymbol{\lambda}}\\[0.4em]

\frac{\partial \boldsymbol{x}_2}{\partial \boldsymbol{\lambda}}
=\boldsymbol{\Phi}(t_2 - t_1) \frac{\partial \boldsymbol{x}_1}{\partial \boldsymbol{\lambda}} +
\boldsymbol{\Gamma}(t_2)\frac{\partial \boldsymbol{I}_2}{\partial \boldsymbol{\lambda}}\\[0.4em]

\frac{\partial \boldsymbol{x}_f}{\partial \boldsymbol{\lambda}} =
\boldsymbol{\Phi}(t_f - t_2) \frac{\partial \boldsymbol{x}_2} {\partial \boldsymbol{\lambda}} +
\boldsymbol{\Gamma}(t_f)\frac{\partial \boldsymbol{I}_f}{\partial \boldsymbol{\lambda}}
$$
Here $I_1$, $I_2$ and $I_f$ correspond to the thrust integral for time $t_0$ to $t_1$, $t_1$ to $t_2$ and $t_2$ to $t_f$, respectively.
It can be seen that the key part of the jacobians are the jacobians of the thrust integral.

As we use the numerical method to evaluate the thrust integral, it can be written as
$$\boldsymbol{i}(\tau)=\begin{pmatrix}
{\boldsymbol{1}_p}_v\cos(\tau)\frac{T}{m(\tau)g_0}\\
{\boldsymbol{1}_p}_v\sin(\tau)\frac{T}{m(\tau)g_0}
\end{pmatrix}$$

$$\begin{aligned}
\boldsymbol{I} &= \int_{t_1}^{t_2} \boldsymbol{i}(\tau)\mathrm{d}\tau\\
&=\frac{t_2-t_1}{N}\sum_{j=0}^4 a_j \boldsymbol{i}(t_1+ j\delta)
\end{aligned}$$
Where ${\boldsymbol{1}_p}_v=\frac{p_v(\tau)}{\|p_v(\tau)\|}$ is the unit vector of $p_v$, and we have
$$\delta = \frac{t_2-t_1}{4}, N=90, a_0=7, a_1=32, a_2=12, a_3=32, a_4=7$$
for Milne’s rule.</br>
Taking derivatives of this expression is quite simple, remember that
$$\boldsymbol{\lambda}(t)=\boldsymbol{\Phi}(t-t_0)\boldsymbol{\lambda}(t_0)$$
The Jacobian of $I$ can be written as
$$
\frac{\partial \boldsymbol{I}}{\partial \boldsymbol{\lambda}}=
\frac{t_2-t_1}{N}\sum_{j=0}^4 a_j
\frac{\partial \boldsymbol{i}(t_1+ j\delta)}{\partial \boldsymbol{\lambda}}
$$
And the Jacobian of $i$ can be written as:
$$
\frac{\partial \boldsymbol{i}(t_1+ j\delta)}{\partial \boldsymbol{\lambda}} =
\frac{T}{m(\tau)g_0}
\begin{pmatrix}
\cos(t_1+ j\delta)I_{3\times3}\\
\sin(t_1+ j\delta)I_{3\times3}
\end{pmatrix}
\frac{\partial {\boldsymbol{1}_p}_v}
    {\partial \boldsymbol{\lambda}}
$$
Now the key part is to calculate the Jacobian
$\frac{\partial {\boldsymbol{1}_p}_v}{\partial \boldsymbol{\lambda}}$.
Apprently, ${\boldsymbol{1}_p}_v$ is  a function of $\boldsymbol{p}_v$, and $\boldsymbol{p}_v$ is a function of $\boldsymbol{\lambda}$.
We can write the Jacobians of the two fuctions:
$$\begin{aligned}
\frac{\partial {{\boldsymbol{1}_p}_v}_i}
    {{\partial p_v}_j}
&=\frac{\|\boldsymbol{p}_v\|-
\frac{{p_v}_i{p_v}_j}{\|p_v\|}}
    {\|\boldsymbol{p}_v\|^2}\\

\frac{\partial {{\boldsymbol{1}_p}_v}}
    {{\partial \boldsymbol{p}_v}}
&=\frac{1}
    {\|\boldsymbol{p}_v\|}
\left(
    1-{\boldsymbol{1}_p}_v{\boldsymbol{1}_p}_v^T
\right)
\equiv \boldsymbol{K}
\end{aligned}
$$
and
$$
\begin{aligned}
    \frac{\partial \boldsymbol{p}_v}
        {\partial \boldsymbol{\lambda}}
    &=
    \begin{pmatrix}
        \frac{\partial \boldsymbol{p}_v}
        {\partial {\boldsymbol{p}_r}_0} &
        \frac{\partial \boldsymbol{p}_v}
        {\partial {\boldsymbol{p}_v}_0}
    \end{pmatrix}\\
    &=
        \begin{pmatrix}
            -\sin(\delta)\boldsymbol{I}_{3\times3} &
            \cos(\delta)\boldsymbol{I}_{3\times3}
        \end{pmatrix}\equiv \boldsymbol{\Lambda}(\delta)\\

    \frac{\partial \boldsymbol{p}_r}{\partial \boldsymbol{\lambda}} &=
    \begin{pmatrix}
        \frac{\partial \boldsymbol{p}_r}
        {\partial {\boldsymbol{p}_r}_0} &
        \frac{\partial \boldsymbol{p}_r}
        {\partial {\boldsymbol{p}_v}_0}
    \end{pmatrix}\\
    &=  
    \begin{pmatrix}
        \cos(\delta)\boldsymbol{I}_{3\times3} &
        \sin(\delta)\boldsymbol{I}_{3\times3}
    \end{pmatrix}
\end{aligned}
$$
where $\delta$ is the time from initial time. Note that the integral lower bound may be different from the initial time as we have multipe stages.
We also give the jacobian of $p_r$ for convenience.

Hence, the jacobian of $I$ can finally be written as

$$
\frac{\partial \boldsymbol{I}}
{\partial \boldsymbol{\lambda}}
=
\frac{t_2-t_1}
    {N}
\sum_{j=0}^4 a_j
\frac{T}
{m(\tau)g_0}
\begin{pmatrix}
    \cos(t_1+ j\delta)\boldsymbol{I}_{3\times3}\\
    \sin(t_1+ j\delta)\boldsymbol{I}_{3\times3}
\end{pmatrix}
\boldsymbol{K}(t_1+ j\delta)\boldsymbol{\Lambda}(j\delta)
$$

#### **1.5.3 State's derivatives respect to $t_f$**

Apart from the Jacobians of the thrust integrals, the state's derivatives respect to $t_f$ are also need to computed the final jacobians. Now we will have a look of the state's expressions.
It can be easily seen that only the third part involved $t_f$, so we can write
$$
\frac{\partial \boldsymbol{x}_f}
    {\partial t_f}=
\frac{\partial \boldsymbol{\Phi}(t_f- t_2)}
    {\partial t_f}x_2+
\frac{\partial \boldsymbol{\Gamma}(t_f)}
    {\partial t_f}\boldsymbol{I}(t_f-t_2)+
\boldsymbol{\Gamma}(t_f)
\frac{\partial \boldsymbol{I}(t_f, t_2)}
    {\partial t_f}
$$
and each part is easy to calculate:
$$
\frac{\partial \boldsymbol{\Phi}(t_f- t_2)}
    {\partial t_f}=
\begin{pmatrix}
    -\sin(t_f- t_2)\boldsymbol{I}_{3\times3} & \cos(t_f- t_2)\boldsymbol{I}_{3\times3}\\
    -\cos(t_f- t_2)\boldsymbol{I}_{3\times3} & -\sin(t_f- t_2)\boldsymbol{I}_{3\times3}
\end{pmatrix}
= -\boldsymbol{\Gamma}(t_f- t_2)\\[0.4em]

\frac{\partial \boldsymbol{\Gamma}(t_f)}
    {\partial t_f}=
\begin{pmatrix}
    \cos(t_f)\boldsymbol{I}_{3\times3} & \sin(t_f)\boldsymbol{I}_{3\times3}\\
    -\sin(t_f)\boldsymbol{I}_{3\times3} & \cos(t_f)\boldsymbol{I}_{3\times3}
\end{pmatrix}
= \boldsymbol{\Phi}(t_f)\\[0.4em]

\frac{\partial \boldsymbol{I}(t_f,t_2)}
    {\partial t_f}=\boldsymbol{i}(t_f)
$$
Thus, the state's derivatives respect to $t_f$ are obtained.

#### **1.5.4 Jacobians of target functions**

After we obtained the Jacobians of the , the jacobian of the target function can be easily calculated.
We write target function $s$ as a function of $\lambda$ and $t_f$, the Jacobian can be written as
$$
\frac{\partial \boldsymbol{s}}
    {\partial \boldsymbol{z}}=
\begin{pmatrix}
    \frac{\partial s_1}
        {\partial \boldsymbol{\lambda}}&
    \frac{\partial s_1}
        {\partial t_f}\\[0.5em]

    \frac{\partial s_2}
        {\partial \boldsymbol{\lambda}}&
    \frac{\partial s_2}
        {\partial t_f}\\

    \vdots & \vdots\\

    \frac{\partial s_7}
    {\partial \boldsymbol{\lambda}}&
    \frac{\partial s_7}
    {\partial t_f}
\end{pmatrix}
$$
Now we will carefully analyze the 7 equations.

let
$$
\begin{pmatrix}
    \frac{\partial \boldsymbol{r}_f}
        {\partial \bullet}\\[0.5em]
    \frac{\partial \boldsymbol{v}_f}
        {\partial \bullet}
\end{pmatrix}=
\frac{\partial \boldsymbol{x}_f}
    {\partial \bullet}
$$
For $s_1$, we have

$$\begin{aligned}
    \frac{\partial s_1}
        {\partial \boldsymbol{\lambda}} &=
    2({r_f}_1 \frac{\partial {r_f}_1}
                   {\partial \boldsymbol{\lambda}}+
    {r_f}_2 \frac{\partial {r_f}_2}
                {\partial \boldsymbol{\lambda}}+
    {r_f}_2 \frac{\partial {r_f}_2}
                 {\partial \boldsymbol{\lambda}})\\
    &=2r_f^T \frac{\partial \boldsymbol{r}_f}
                  {\partial \boldsymbol{\lambda}}\\
    \frac{\partial s_1}
    {\partial t_f}&=
    2r_f^T \frac{\partial r_f}
                {\partial t_f}
\end{aligned}
$$
For $s_2$, $s_3$ and $s_4$, they only contain ${v_f}_1$, ${v_f}_2$ and ${v_f}_3$, respectively.We can thus acquire the following results:

$$
\begin{pmatrix}
    \displaystyle \frac{\partial s_2}
                    {\partial \boldsymbol{\lambda}}&
    \displaystyle \frac{\partial s_2}
                    {\partial t_f}\\[1.3em]

    \displaystyle \frac{\partial s_3}
                    {\partial \boldsymbol{\lambda}}&
    \displaystyle \frac{\partial s_3}
                    {\partial t_f}\\[1.3em]

    \displaystyle \frac{\partial s_4}
                    {\partial \boldsymbol{\lambda}}&
    \displaystyle \frac{\partial s_4}
                    {\partial t_f}
\end{pmatrix} = \boldsymbol{I}_{3\times3}
\begin{pmatrix}
    \displaystyle \frac{\partial v_f}
                    {\partial \boldsymbol{\lambda}} &
    \displaystyle \frac{\partial v_f}
                    {\partial t_f}
\end{pmatrix}
$$
The last three equations' derivate can be obtained with similar process, here we gave the results without detailed analysis.
$$
\frac{\partial s_5}
    {\partial \boldsymbol{\lambda}}=
v_f^T\frac{\partial {\boldsymbol{p}_r}_f}
        {\partial \boldsymbol{\lambda}}+
{p_r}_f^T\frac{\partial \boldsymbol{v}_f}
        {\partial \boldsymbol{\lambda}}-
r_f^T\frac{\partial {\boldsymbol{p}_v}_f}
        {\partial \boldsymbol{\lambda}}-
{p_v}_f^T\frac{\partial \boldsymbol{r}_f}
            {\partial \boldsymbol{\lambda}}+
\frac{T_{max}}
    {m_f g_0}
\frac{{\boldsymbol{p}_v}_f^T}
    {\|{\boldsymbol{p}_v}_f\|}
\frac{\partial {\boldsymbol{p}_v}_f}
{\partial \boldsymbol{\lambda}}\\[0.4em]

\frac{\partial s_5}
    {\partial t_f}=
v_f^T \frac{\partial {\boldsymbol{p}_r}_f}
        {\partial t_f}+
{p_r}_f^T\frac{\partial \boldsymbol{v}_f}
            {\partial t_f}-
r_f^T \frac{\partial {\boldsymbol{p}_v}_f}
            {\partial t_f}-
{p_v}_f^T\frac{\partial \boldsymbol{r}_f}
            {\partial t_f}+
\frac{T_{max}}
    {m_f g_0}
\frac{{\boldsymbol{p}_v}_f^T}
    {\|\boldsymbol{p}_v\|}
\frac{\partial {\boldsymbol{p}_v}_f}
    {\partial t_f}+
\frac{T_{max}^2}
    {m_f^2 g_0 v_{ex}}
\|{\boldsymbol{p}_v}_f\|\sqrt{\frac{r_0}{g_0}}\\[0.4em]

\frac{\partial s_{6,7}}
    {\partial \boldsymbol{\lambda}}=
\left\{
    [{\boldsymbol{p}_r}_f+2\kappa(\boldsymbol{r}_f-\boldsymbol{r}^*)]^T \boldsymbol{A}_i
    \frac{\partial \boldsymbol{r}_f}
        {\partial \boldsymbol{\lambda}}+\boldsymbol{r}_f^T \boldsymbol{A}_i^T
    (
        \frac{\partial {\boldsymbol{p}_r}_f}
            {\partial \boldsymbol{\lambda}}+
        2\kappa\frac{\partial \boldsymbol{r}_f}
                    {\partial \boldsymbol{\lambda}}
    )
\right\}\\[0.4em]

\frac{\partial s_{6,7}}
{\partial t_f}=
\left\{
    [{\boldsymbol{p}_r}_f+2\kappa(\boldsymbol{r}_f-\boldsymbol{r}^*)]^T \boldsymbol{A}_i
    \frac{\partial \boldsymbol{r}_f}
        {\partial t_f}+\boldsymbol{r}_f^T \boldsymbol{A}_i^T
    (
        \frac{\partial {\boldsymbol{p}_r}_f}
            {\partial t_f}+
        2\kappa\frac{\partial \boldsymbol{r}_f}
                    {\partial t_f}
    )
\right\}
$$
where $\boldsymbol{A} \boldsymbol{r}_f = \boldsymbol{y}$

Now the root finding problem is ready to be solved.

## 2.AAPDG

The augmented APDG algorithm is capable to perform landing without suffering solved fail to converge problems.

### 2.1 E-guidance and APDG

If we assume the gravity is contant, then the state equations will be

$$
\dot{\boldsymbol{r}} = \boldsymbol{v}\\
\dot{\boldsymbol{v}} = \boldsymbol{g} + \boldsymbol{a_T}\\
$$
we want to arrive the target postion with target speed in the given time, that is,

$$
\boldsymbol{r}(t_f) = \boldsymbol{r}^*\\
\boldsymbol{v}(t_f) = \boldsymbol{v}^*\\
$$
In E-guidance, we assume the acceleration is linear with time,

$$
\boldsymbol{a}_T=\boldsymbol{c_0}+\boldsymbol{c_1}t
$$
coefficients $\boldsymbol{c_0}$ and $\boldsymbol{c_1}$ can be determined through satisfying the terminal conditions, and we'll finally arrive at the E-guidance law:

$$
\boldsymbol{a}_T=
-\frac{2}
    {t_{go}}
[\boldsymbol{v}^*-\boldsymbol{v}]+
\frac{6}{t_{go}^2}[\boldsymbol{r}^*-\boldsymbol{r}-\boldsymbol{v}t_{go}]-\boldsymbol{g}
$$
where $t_{go} = t_f-t$ is the time till landing.
If we use a quadratic acceleration approximation, with the final acceleration condition
$\boldsymbol{a}_T(t_f)=\boldsymbol{a}_T^*$, one can derive the APDG law:

$$
\boldsymbol{a}_T=
-\frac{6}
    {t_{go}}
[\boldsymbol{v}^*-\boldsymbol{v}]+
\frac{12}{t_{go}^2}[\boldsymbol{r}^*-\boldsymbol{r}-\boldsymbol{v}t_{go}]+\boldsymbol{a}_T^*
$$
It is known that E-guidance is more propellant optimal than APDG, but the latter gives ability to conform the final acceleration condition.

### 2.2 AAPDG law

It can be observed that E-guidance and APDG have a similar structure, with different coefficients. If we assume a implicit guidance law

$$
\boldsymbol{a}_T= \boldsymbol{a}_d
-\frac{k_v}
    {t_{go}}
[\boldsymbol{v}^*-\boldsymbol{v}]-
\frac{k_r}{t_{go}^2}[\boldsymbol{r}^*-\boldsymbol{r}]
$$
and assume the acceleration to be linear, the coefficient of jerk can be eliminate by setting $k_v=\frac{1}{3}k_r+2$, and gives following guidance law:

$$
\boldsymbol{a}_T=
\frac{2}
    {t_{go}}
(1-\frac{k_r}{3})
[\boldsymbol{v}^*-\boldsymbol{v}]+
\frac{k_r}{t_{go}^2}[\boldsymbol{r}^*-\boldsymbol{r}-\boldsymbol{v}t_{go}]+\frac{1}{6}(k_r-6)\boldsymbol{a}_T^*+\frac{1}{6}(k_r-12)\boldsymbol{g}
$$
when $k_r=6$, this gives E-guidance law, and when $k_r=12$, this gives APDG law. Any value between 6 and 12 is feasible, and larger $k_r$ typically uses more propellant, but gives a better trajectory.

However, the final acceleration will be different from $\boldsymbol{a}_T^*$ except the $k_r=12$ case. When close enough to the target, the guidance law gives

$$
\boldsymbol{a}_T(t_f)=(\frac{k_r}{6}-1)\boldsymbol{a}_T^*+(\frac{k_r}{6}-2)\boldsymbol{g}
$$
typically we want the final acceleration be $\boldsymbol{a}_T^*=-\kappa\boldsymbol{g}$, where $\kappa > 1$. Consequently, the final acceleration will be
$$
\boldsymbol{a}_T(t_f)=-\left[(\kappa-1)(\frac{k_r}{6}-1)+1\right]\boldsymbol{g}
$$
and we can calculate the correct $\kappa$ value to make it conform with our demand.

For further discussion see the original paper for details.
