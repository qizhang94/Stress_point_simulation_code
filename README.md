# Onsite talk subtitle

## THM 1
For this research topic, we have made novel contributions to the coupling theory. 

In terms of the mass balance equation, which is also known as the fluid seepage equation, we start from its primitive form shown in the first row. By using the material time derivative, the Darcy velocity $q_f$ comes into play. To get rid of the porosity evolution term ${\rm d} \phi/ {\rm d} t$, we adopt the principle of superposition, axiomatic relation between partial and total stress tensors, and material state equation, so we derive an equation for ${\rm d} \phi/ {\rm d} t$. The combination of the previous two equations finally give our novel result. 

Here, the novelty lies in the new general volumetric thermal expansion coefficient beta e, which incorporates material anisotropy. In addition, under the case of isotropy, this beta e is exactly the same as that of Cheng (2016), who derived this coefficient using a different method. Besides, the coefficient $S$ is also the most general storage coefficient in poromechanics.

## THM 2
In terms of the thermal equation, which is also known as the energy balance equation, we start from the balance of energy equation for each constituent. It says that the material time derivative of the internal energy is equal to the summation of power term, energy source, energy exchange between solid and fluid, minus the time derivative of the kinetic energy.

After a series of simplification, we finally obtain the novel balance of energy equation for the whole mixture, which is shown here. Compare with the uncoupled standard thermal equation, our equation contains several new energy source terms, as shown here.

How to quantify these new energy source terms? Well, the scaling analysis is a powerful technique! The results are summarized here. It turns out that all the new terms have the same order of magnitude. By taking the ratio, we produce a new dimensionless number omega_T and a new rule of thumb! In this dimensionless number, k_f is the characteristic fluid mobility, k_T is the characteristic thermal conductivity, Delta P is the characteristic fluid pressure variation, and Delta T is the characteristic temperature variation. If this omega_T is way less than 1, then for the thermal equation, we could use the uncoupled form, despite the coupling exists.

## THM 3
For the numerical scheme, we inherit from the so-called node-based smoothed finite element method, or NS-FEM for short. Compare with the standard FEM, the NS-FEM assembles the system residual and tangent operator on the smoothing domain, as a result, it avoids the variable mapping when there is re-mesh, because all the information is stored at the nodes.

As you can see from the left figure, after re-mesh, we need to map the stress, strain, and other history variables from the old Gauss points to new Gauss points, which will lead to some errors. While in NS-FEM, the nodes are preserved, and all the information moves with the node. No accuracy is lost during this process.

Nevertheless, the standard NS-FEM is not stable enough, and therefore, we adopt the gradient-based stabilization method to recover the strain gradient within the smoothing domain. In this way, the method becomes more robust.

## THM 4
Here, we show one verification example, which is the thermo-elastic consolidation. In this example, we apply a surface load and a higher temperature on top of the soil column. As a result, the soil column would have a gradual settlement, followed by a surface rebound. 

It is worth to mention that in our setting, the surface load and higher temperature are not applied through the boundary condition, instead, we apply through the multiphysics contact. We assume there is a drained rigid indenter on top of the soil column, and it is in steady-state for pore pressure and temperature fields. The surface load and temperature are applied on the indenter, which would then transfer to the soil column. In that way, we also verify the implementation of the contact algorithm.

A good agreement could be observed with the reference results.

## THM 5
We apply our method to investigate the pipeline penetration process considering thermal effects. The geometry and mesh of the submarine pipeline and seabed soil is shown here on the left.

In this figure, we compare the simulation results with experimental findings regarding soil configuration during pipeline penetration. For stiff soil with a large Young’s modulus, soil heave occurs after pipeline penetration, consistent with the findings of Dingle et al. Conversely, for soft clay with a small Young’s modulus, no soil heave is observed, and the pipeline sinks into the soil, aligning with the results of the model test.

On the right, we show the evolution of equivalent plastic strain and excess pore pressure during the pipeline penetration. On the bottom, a video depicts the adaptive remesh process, we can observe that the mesh quality remains satisfactory during pipeline penetration.

## THM 6
Now we look at some quantitative results on the penetration resistance. The left column suggests that penetration resistance is insensitive to temperature for frictionless contact. However, for frictional contact, as shown by the middle column, a higher pipeline temperature leads to a higher excess pore water pressure. The presence of excess pore water pressure reduces the effective normal force on the soil, thereby lowering friction, which explains why the penetration resistance of the 95℃ pipeline is relatively lower. 

Additionally, the increase in the thermal expansion coefficient enhances the expansion tendency of the soil, consequently increasing resistance, as shown by the green line on the right column, compare with that in the middle.

Therefore, the thermal effects on resistance can be summarized as a competition mechanism between friction degradation and thermal expansion.

## CLOSURE 1
First of all, I like to highlight the availability of open-source computer codes specifically designed for point simulation of shale and the SNS-PFEM method. These codes are accessible to the public, promoting transparency and collaboration within the research community.

The repository is well-documented and actively maintained, making it a valuable resource for researchers in this field.

## CLOSURE 2
There are some main takeaways of this talk.

