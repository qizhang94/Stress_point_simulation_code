## Some important remarks about the MCC_UMAT.m file
- This file cannot be perfect, it would lead to non-convergence issue for some (in fact, **a lot of** 😄) corner cases, however, for basic triaxial test, it is enough, and the typical results are given [here](https://mp.weixin.qq.com/s/HGIRHsGQu3TsCcgKJAlIFg). As a comparison, the triaxial test using Drucker-Prager model could be visualized [here](https://github.com/qizhang94/Stress_point_simulation_code/blob/main/Drucker-Prager/Drucker-Prager_triaxial.png).
- Some possible improvements:
    - Try to change `sigma_iter = sigma_old` to `sigma_iter = sigma_trial`
    - Use `while norm(R)/r0 >= error_tol`, _i.e._, use the residual of totally 11 equations, instead of the first 9 equations
    - Usually, the reason for non-convergence is from the inversion of the _J_ matrix of the return-mapping algorithm, so one may need to check `rcond(J)` in MATLAB to ensure it is not too small
    - Also, because we have an exponential update for the $p_c$ (preconsolidation pressure/stress), so if the plastic multiplier $\Delta \lambda$ is a large positive or large negative number, that would be a disaster. Therefore, please check `abs(Dlambda_iter)` in MATLAB
- Below is another non-convergence issue, but the reason is just my own guess 😂
    - At the beginning of a time step, when the strain increment $\Delta \epsilon$ is zero, however, it is surprisingly to find that the statement `while norm(R)/r0 >= error_tol` is true!! (because theoretically, at the beginning of a time step, this statement should be false!!) As a result, the local Newton iteration would never converge given the strain increment is zero
    - We think it is because the `sigma_iter` got another update or change after the convergence (`norm(R)/r0 < error_tol`) in the **previous time step**. This update could make `norm(R)/r0 >= error_tol` in the current step, although the violation is very very small (**This is only our guess**)
