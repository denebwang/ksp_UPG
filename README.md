# ksp_UPG

<<<<<<< HEAD
universal powered guidance in Kerbal Space Program based on krpc

## 1. About UPG

Universal Powered Guidance is a versatile powered guidance developed by Lu, et al. [1-3] , which is based on optimal control theory.

By setting different constraints, the algorithms can be used for differnet purpose including landing, ascent, etc. Here only the powered descent is implemented.


=======
currently, only powered descent module is implemented.
## How to use
In `target.json`, enter the target parameters.
Especially, if you want just to land somewhere rather than to target, set `k` to 0.
Run `Bolza-Landing.py` to get the guidance. The program will handle everything.
>>>>>>> 9ef79e3d5a96d53379fa8f3bc95419f1f6ae55bb

## References

[1] Lu, Ping. "Propellant-optimal powered descent guidance." Journal of Guidance, Control, and Dynamics 41.4 (2018): 813-826.

[2] Lu, Ping, et al. "Rapid optimal multiburn ascent planning and guidance." Journal of Guidance, Control, and Dynamics 31.6 (2008): 1656-1664.

[3] Lu, Ping, Stephen Forbes, and Morgan Baldwin. "A versatile powered guidance algorithm." AIAA Guidance, Navigation, and Control Conference. 2012.
