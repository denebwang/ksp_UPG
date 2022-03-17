# ksp_UPG

universal powered guidance in Kerbal Space Program based on krpc

## About UPG

Universal Powered Guidance is a versatile powered guidance developed by Lu, et al. [1-3] , which is based on optimal control theory.

By setting different constraints, the algorithms can be used for differnet purpose including landing, ascent, etc. Here only the powered descent is implemented.

## Using UPG

It is preferred to use the `Bolza-Landing.py` file to perform landing guidance.
To begin with, specify the target latitude and longitude in the parameter section, and you may need to adjust some of the parameters to get it working properly.

UPG descent module works bad from orbit, you may want to get in a descenting trajectory ( e.g. pe=15km, above target) when landing for moon, and for mars you can use upg after the airbreak phase (this is what upg is originally developed for).

There is also a file to perform landing purely by the Augmented Apollo Powered Descent Guidance (AAPDG) [4]. This algorithm is used for UPG's terminal guidance, and itself can be used to perform  powered descent guidance. AAPDG doesn't suffer from the numerically difficult problem, which cause the UPG fail to converge when getting closer to the target.

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