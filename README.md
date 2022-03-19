# ksp_UPG

universal powered guidance in Kerbal Space Program based on krpc. The solver doesn't work very well, for more stable performance you can use the Apollo guidance.

## About UPG

Universal Powered Guidance is a versatile powered guidance developed by Lu, et al. [1-3] , which is based on optimal control theory.

By setting different constraints, the algorithms can be used for differnet purpose including landing, ascent, etc. Here only the powered descent is implemented.

For the theory, refer to the paper list at the end. ~~There is also a summary in `theory.ipynb`. Notebook is used as github cannot show formulae.~~ I can't get github preview working right, if you want to see the theories you can download the [theory.md](https://github.com/denebwang/ksp_UPG/blob/master/theory.md) file.

## Using UPG

It is preferred to use the `Bolza-Landing.py` file to perform landing guidance.
To begin with, specify the target latitude and longitude in the parameter section, and you may need to adjust some of the parameters to get it working properly.

### Step 1: Preparation

UPG descent module works bad from orbit, you may want to get in a descenting trajectory (e.g. pe=15km, above target) when landing for moon, and for mars you can use upg after the airbreak phase (this is what upg is originally developed for).

Besides, UPG does not support multi-stage landers. You can burn out the previous to get into a landing trajectory, and then use the UPG for fine tuning.

Once you have entered the trajectory, you can initialize the UPG. Modify the parameters, especially the `lon` and `lat` ones. Other parameters may vary depending on craft and the body. The UPG will automatically estimate the optimal "coast" time. Here "coast" does not mean shutdown the engine, it will put throttle to minimum to avoid ignition.

### Step 2: Wait for start

UPG comes with powered descent initial(PDI) phase, that is, decide when to start the descent. This is done by running the softlanding method, and compare the downrange of softlanding and actual downrange. The descent will start when it slightly overshoots.

### Step 3: Main guidance

After the main algorithm kicks in, the craft will be guided towards the target. You may notice that the solver will fail to converge at times, and especially when close to target. This is normal, and should not be a problem as long as you can finally get there. The guidance will use the previously converged solution for open loop guidance.

### Step 4: Terminal guidance

As mentioned above, the solver will not converge well near the target. As a result, a terminal guidance of Augmented Apollo Powered Descent Guidance (AAPDG) [4] is added. There is also a file to perform landing purely by AAPDG, but use UPG to decide start time and total time for landing. AAPDG doesn't suffer from the numerically difficult problem, which cause the UPG fail to converge when getting closer to the target. AAPDG itself enconters singularity when time to-go approches 0, so a fixed speed descent is added at last.

## What about efficiency?

Here's a comparison of the actual delta-v comsumption of these 2 algorithms, when performing lunar landing from a 100km@15km oribit. Note that the delta-v of upg is always greater than what it says, as the terminal guidance will use more.

|algorithm   |delta-v|
|------------|-------|
|UPG         |1790   |
|AAPDG, kr=6 |1835   |
|AAPDG, kr=8 |1847   |
|AAPDG, kr=10|1866   |
|AAPDG, kr=12|1913   |

The kr = 6 case is equivalent to E-guidance, and the kr = 12 case is equivalent to original APDG used for Apollo landers. It is not guaranteed that the both algorithm will not crash into terrain. Nor can AAPDG limit the thrust within the craft's capabilities, it can be saturated at times. Sometimes this can lead to crashed into terrain.

## Dependencies

numpy  
scipy  
krpc  
colorama (for display)

## References

[1] P. Lu, “Propellant-optimal powered descent guidance,” J. Guid. Control Dyn., vol. 41, no. 4, pp. 813–826, 2018, doi: 10.2514/1.G003243.

[2] P. Lu, B. J. Griffin, G. A. Dukeman, and F. R. Chavez, “Rapid optimal multiburn ascent planning and guidance,” in Journal of Guidance, Control, and Dynamics, 2008, vol. 31, no. 6, pp. 1656–1664. doi: 10.2514/1.36084.

[3] P. Lu, S. Forbes, and M. Baldwin, “A versatile powered guidance algorithm,” 2012. doi: 10.2514/6.2012-4843.

[4] P. Lu, “Augmented Apollo Powered Descent Guidance,” J. Guid. Control Dyn., vol. 42, no. 3, pp. 447–457, 2019, doi: 10.2514/1.G004048.
