# Theory

## 2.1 State equations

### 2.2 Constrains for root finding problem

### 2.3 Jacobians

In this part, the Jacobians of the target function are given. Note that the constraints can be analytically calculated, the jacobians can also be derived.

#### 2.3.1 Summary of the equations

For convenience, we will breifly summarize how the target function computed.\

Given the target state $\boldsymbol{r}^*$ and $\boldsymbol{v}^*$, we should first calculate the final state $\boldsymbol{x}$ corresponds to the unknowns $\boldsymbol{z} = (\boldsymbol{p}_v^T \; \boldsymbol{p}_r^T \; t_f)^T \in R^7$. As the thrust will be 3-arc, we will have to calculate the states one by one.

Firstly, we calculate the state at time $t_1$:
$$
\boldsymbol{x}_1 =
\boldsymbol{\Phi}(t_1 - t_0) \, x_0 +
\boldsymbol{\Gamma}(t_1) \, \boldsymbol{I}(t_1 - t_0)
$$

Then we do the same thing for $t_2$ and $t_f$:

$$
\boldsymbol{x}_2 =
\boldsymbol{\Phi}(t_2 - t_1) \, x_1 +
\boldsymbol{\Gamma}(t_2) \, \boldsymbol{I}(t_2 - t_1)\\

\boldsymbol{x}_f =
\boldsymbol{\Phi}(t_f - t_2) \, x_2 +
\boldsymbol{\Gamma}(t_f) \, \boldsymbol{I}(t_f - t_2)
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

#### **2.3.2 Jacobians of thrust integral**

Now that we have obtained the target function $s=(s_1, s_2, s_3, s_4, s_5, s_6, s_7)^T$, we have to calculate the Jacobian
$$
\boldsymbol{J}=\frac{\partial s}{\partial z}
$$
To simplify the results, we use
$$\lambda = \begin{pmatrix}
p_r\\
P_v
\end{pmatrix}$$
to replace $p_v$ and $p_r$ in  $z$.
If we derivate $x$ against $\lambda$, we can write the results one by one:
$$\frac{\partial x_1}{\partial \lambda} = \Gamma(t_1)
\frac{\partial I_1}{\partial \lambda}$$

$$\frac{\partial x_2}{\partial \lambda} = \Phi(t_2 - t_1)\frac{\partial x_1}{\partial \lambda} +
\Gamma(t_2)\frac{\partial I_2}{\partial \lambda}$$

$$\frac{\partial x_f}{\partial \lambda} = \Phi(t_f - t_2)\frac{\partial x_2}{\partial \lambda} +
\Gamma(t_f)\frac{\partial I_f}{\partial \lambda}$$
Here $I_1$, $I_2$ and $I_f$ correspond to the thrust integral for time $t_0$ to $t_1$, $t_1$ to $t_2$ and $t_2$ to $t_f$, respectively.
It can be seen that the key part of the jacobians are the jacobians of the thrust integral.</br>
As we use the numerical method to evaluate the thrust integral, and can be written as
$$i(\tau)=\begin{pmatrix}
{\boldsymbol{1}_p}_v\cos(\tau)\frac{T}{m(\tau)g_0}\\
{\boldsymbol{1}_p}_v\sin(\tau)\frac{T}{m(\tau)g_0}
\end{pmatrix}$$

$$\begin{aligned}
I &= \int_{t_1}^{t_2} i(\tau)\mathrm{d}\tau\\
&=\frac{t_2-t_1}{N}\sum_{j=0}^4 a_j i(t_1+ j\delta)
\end{aligned}$$
Where ${\boldsymbol{1}_p}_v=\frac{p_v(\tau)}{\|p_v(\tau)\|}$ is the unit vector of $p_v$, and we have
$$\delta = \frac{t_2-t_1}{4}, N=90, a_0=7, a_1=32, a_2=12, a_3=32, a_4=7$$
for Milne’s rule.</br>
Taking derivatives of this expression is quite simple, remember that
$$\lambda(t)=\Phi(t-t_0)\lambda(t_0)$$
The Jacobian of $I$ can be written as
$$
\frac{\partial I}{\partial \lambda}=\frac{t_2-t_1}{N}\sum_{j=0}^4 a_j
\frac{\partial i(t_1+ j\delta)}{\partial \lambda}
$$
And the Jacobian of $i$ can be written as:
$$
\frac{\partial i(t_1+ j\delta)}{\partial \lambda} = \frac{T}{m(\tau)g_0}
\begin{pmatrix}
\cos(t_1+ j\delta)I_{3\times3}\\
\sin(t_1+ j\delta)I_{3\times3}
\end{pmatrix}
\frac{\partial {\boldsymbol{1}_p}_v}{\partial \lambda}
$$
Now the key part is to calculate the Jacobian
$\frac{\partial {\boldsymbol{1}_p}_v}{\partial \lambda}$.
Apprently, ${\boldsymbol{1}_p}_v$ is  a function of $p_v$, and $p_v$ is a function of $\lambda$.
We can right the Jacobians of the two fuctions:
$$\begin{aligned}
\frac{\partial {{\boldsymbol{1}_p}_v}_i}{{\partial p_v}_j}
&=\frac{\|p_v\|-\frac{{p_v}_i{p_v}_j}{\|p_v\|}}{\|p_v\|^2}\\
\frac{\partial {{\boldsymbol{1}_p}_v}}{{\partial p_v}}
&=\frac{1}{\|p_v\|}\left(1-{\boldsymbol{1}_p}_v{\boldsymbol{1}_p}_v^T\right)
\equiv K
\end{aligned}
$$
and
$$
\begin{aligned}
\frac{\partial p_v}{\partial \lambda}&=
    \begin{pmatrix}
    \frac{\partial p_v}{\partial {p_r}_0} & \frac{\partial p_v}{\partial {p_v}_0}
    \end{pmatrix}\\
&=\begin{pmatrix}
    -\sin(\delta)I_{3\times3} & \cos(\delta)I_{3\times3}
    \end{pmatrix}\equiv \Lambda(\delta)\\
\frac{\partial p_r}{\partial \lambda}&=
    \begin{pmatrix}
    \frac{\partial p_r}{\partial {p_r}_0} & \frac{\partial p_r}{\partial {p_v}_0}
    \end{pmatrix}\\
&=\begin{pmatrix}
    \cos(\delta)I_{3\times3} & \sin(\delta)I_{3\times3}
    \end{pmatrix}
\end{aligned}
$$
where $\delta$ is the time from initial time. Note that the integral lower bound may be different from the initial time as we have multipe stages.
We also give the jacobian of $p_r$ for convenience.</br>
Hence, the jacobian of $I$ can finally be written as
$$
\frac{\partial I}{\partial \lambda} =\frac{t_2-t_1}{N}\sum_{j=0}^4 a_j
\frac{T}{m(\tau)g_0}
    \begin{pmatrix}
    \cos(t_1+ j\delta)I_{3\times3}\\
    \sin(t_1+ j\delta)I_{3\times3}
    \end{pmatrix}
K(t_1+ j\delta)\Lambda(j\delta)
$$

#### **2.3.3 State's derivatives respect to $t_f$**

Apart from the Jacobians of the thrust integrals, the state's derivatives respect to $t_f$ are also need to computed the final jacobians. Now we will have a look of the state's expressions.
It can be easily seen that only the third part involved $t_f$, so we can just write
$$
\frac{\partial x_f}{\partial t_f}=\frac{\partial \Phi(t_f- t_2)}{\partial t_f}x_2+\frac{\partial \Gamma(t_f)}{\partial t_f}I(t_f-t_2)+\Gamma(t_f)\frac{\partial I(t_f, t_2)}{\partial t_f}
$$
and each part is easy to calculate:
$$
\frac{\partial \Phi(t_f- t_2)}{\partial t_f}=
\begin{pmatrix}
-\sin(t_f- t_2) & \cos(t_f- t_2)\\
-\cos(t_f- t_2) & -\sin(t_f- t_2)
\end{pmatrix}
= -\Gamma(t_f- t_2)
$$

$$
\frac{\partial \Gamma(t_f)}{\partial t_f}=
\begin{pmatrix}
\cos(t_f) & \sin(t_f)\\
-\sin(t_f) & \cos(t_f)
\end{pmatrix}
= \Phi(t_f)
$$

$$
\frac{\partial I(t_f,t_2)}{\partial t_f}=i(t_f)
$$
Thus, the state's derivatives respect to $t_f$ are obtained.

#### **2.3.4 Jacobians of target functions**

After we obtained the Jacobians of the , the jacobian of the target function can be easily calculated.
We write target function $s$ as a function of $\lambda$ and $t_f$, the Jacobian can be written as
$$
\frac{\partial s}{\partial z} =
\begin{pmatrix}
\frac{\partial s_1}{\partial \lambda} & \frac{\partial s_1}{\partial t_f}\\[0.5em]
\frac{\partial s_2}{\partial \lambda} & \frac{\partial s_2}{\partial t_f}\\
\vdots & \vdots\\
\frac{\partial s_7}{\partial \lambda} & \frac{\partial s_7}{\partial t_f}
\end{pmatrix}
$$
Now we will carefully analyze the 7 equations.</br>
let
$$
\begin{pmatrix}
\frac{\partial r_f}{\partial \bullet}\\[0.5em]
\frac{\partial v_f}{\partial \bullet}
\end{pmatrix}=\frac{\partial x_f}{\partial \bullet}
$$
For $s_1 = r_f^T r_f - {r^*}^Tr^*$, we have
$$\begin{aligned}
\frac{\partial s_1}{\partial \lambda} &= 2({r_f}_1 \frac{\partial {r_f}_1}{\partial \lambda}+
{r_f}_2 \frac{\partial {r_f}_2}{\partial \lambda}+{r_f}_2 \frac{\partial {r_f}_2}{\partial \lambda})\\
&=2r_f^T\frac{\partial r_f}{\partial \lambda}
\end{aligned}
$$
and
$$
\frac{\partial s_1}{\partial t_f} = 2r_f^T\frac{\partial r_f}{\partial t_f}
$$
For $s_2$, $s_3$ and $s_4$, they only contain ${v_f}_1$, ${v_f}_2$ and ${v_f}_3$, respectively.We can thus acquire the following results:
$$
\begin{pmatrix}
\displaystyle\frac{\partial s_2}{\partial \lambda}
& \displaystyle\frac{\partial s_2}{\partial t_f}\\[1.3em]
\displaystyle\frac{\partial s_3}{\partial \lambda}
& \displaystyle\frac{\partial s_3}{\partial t_f}\\[1.3em]
\displaystyle\frac{\partial s_4}{\partial \lambda}
& \displaystyle\frac{\partial s_4}{\partial t_f}
\end{pmatrix} =
I_{3\times3}
\begin{pmatrix}
\displaystyle
\frac{\partial v_f}{\partial \lambda} &
\displaystyle
\frac{\partial v_f}{\partial t_f}
\end{pmatrix}
$$
The last three equations' derivate can be obtained with similar process, here we gave the results without detailed analysis.
$$
\frac{\partial s_5}{\partial \lambda} =
v_f^T\frac{\partial {p_r}_f}{\partial \lambda} + {p_r}_f^T\frac{\partial v_f}{\partial \lambda} - r_f^T\frac{\partial {p_v}_f}{\partial \lambda} - {p_v}_f^T\frac{\partial r_f}{\partial \lambda} +
\frac{T_{max}}{m_f g_0} \frac{{p_v}_f^T}{\|{p_v}_f\|} \frac{\partial {p_v}_f}{\partial \lambda}\\

\frac{\partial s_5}{\partial t_f} =
v_f^T \frac{\partial {p_r}_f}{\partial t_f} + {p_r}_f^T\frac{\partial v_f}{\partial t_f} -
r_f^T \frac{\partial {p_v}_f}{\partial t_f} - {p_v}_f^T\frac{\partial r_f}{\partial t_f}+\frac{T_{max}}{m_f g_0}\frac{{p_v}_f^T}{\|p_v\|}\frac{\partial {p_v}_f}{\partial t_f}+\frac{T_{max}^2}{m_f^2 g_0 v_{ex}}\|{p_v}_f\|\sqrt{\frac{r_0}{g_0}}
$$

$$
\frac{\partial s_{6,7}}{\partial \lambda}=
\left\{[{p_r}_f+2\kappa(r_f-r^*)]^T A_i \,\frac{\partial r_f}{\partial \lambda}+
r_f^T A_i^T (\frac{\partial {p_r}_f}{\partial \lambda}+2\kappa\frac{\partial r_f}{\partial \lambda})\right\}\\

\frac{\partial s_{6,7}}{\partial t_f}=
\left\{[{p_r}_f+2\kappa(r_f-r^*)]^T A_i \,\frac{\partial r_f}{\partial t_f}+
r_f^T A_i^T(\frac{\partial {p_r}_f}{\partial t_f}+2\kappa\frac{\partial r_f}{\partial t_f})\right\}
$$
