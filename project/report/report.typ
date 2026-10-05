#import "template/frontpage.typ": conf
#import "template/components.typ": *

#import "@preview/wrap-it:0.1.1": *

#show: conf.with(
  lang: "en",
  cours: "SF2529 - Inverse Problems",
  subject: "Project",
  title: "Tomographic Reconstruction",
  students: (
    (name: "Lucas Ahou", noma: "Group 8"),
  ),
  teachers: (
    (name: "C. Felipe Jerez Hanckes"),
    (name: "A. Janson"),
    (name: "J. Krook"),
  ),
  heading_numbering: "1.a" // any pattern or none if you don't want to number
)

#let xsol = $x^delta_lambda$

= Introduction

This project aims to use regularization techniques to reconstruct a 3D CT scan of a human head. A CT scan projects X-rays through an object following a given ray pattern (here, a cone beam was used) and a sensor at the other end receives those attenuated rays. This is done for many angles which produces a sinogram, i.e. a graph with the projection angle on the abscissa axis and the values of the incoming rays at the sensor on the ordinate axis. In practice, this sinogram is affected by noise. The reconstruction problems is therefore, given a sinogram data $y$, to recover an image $x$ that is close to the original one.

The problem is that the forward projection operator tends to smooth noise. Intuitively, its inverse will thus have the opposite effect, i.e. amplify noise.
While a technique like least-square approximation will minimize the residual between the given data and the image of our solution, it will try to fit the noise which is not what we want in this application.\
A better solution is thus to use a regularizer to better capture some features of the images that we would like to recover. The mathematical problem we would like to solve is the following one:

$
  xsol := min_x 1/2 ||A x - y^delta|| + lambda H_(epsilon)(gradient x)
$

Here, the regularizer is the Hubert function defined as:

$
  H_(epsilon)(z) := sum_i h_(epsilon)(z_i)
$

with:

$
  h_(epsilon)(t) := 
  cases(
    t^2/(2 epsilon) quad &"if" |t| <= epsilon,
    |t| - epsilon/2 quad &"otherwise"
  )
$

We could have use the total variation (TV) regularizer but it is not differentiable and thus harder to optimize. The Hubert function is a smooth approximation of the TV regularizer that makes use of a parameter $epsilon$ to control the smoothness of the approximation.
With this, we can use a gradient descent method to find the solution of our problem. The parameter $lambda$ is used to control the trade-off between fitting the data and smoothing the noise. A small value of $lambda$ will lead to a solution that fits the data well but is noisy, while a large value of $lambda$ will lead to a smooth solution that does not fit the data well.

In the following sections, I will first solve this problem on the 2D Shepp-Logan phantom for testing purposes and to tune the method before applying it to the 3D CT scan of a human head, and compare with the filtered backprojection (FBP) method.

= Reconstruction in 2D

In this section, I will present the parameter choices and the results of the reconstruction of the 2D Shepp-Logan phantom.

== Step size

As mentioned in the introduction, gradient descent is fitted to solve this problem as the objective is differentiable, convex and smooth. Gradient descent is guaranteed to converge for a convex, L-smooth function if the step size satisfy $0 < t < 2/L$. In our case, the first term of the objective has a smoothness constant $L_1 = ||A||^2$ and the second term has a smoothness constant $L_2 = 1/epsilon$#footnote[Given in assignment]. Composed with the gradient operator and the $lambda$ scalar, this gives:

$
  L = ||A||^2 + lambda (||gradient||^2)/(epsilon)
$


In the code, I chose to use a step size $t = 1.9/L$.

== Choice of $epsilon$ and $lambda$

First, I wanted to determine a good value for $epsilon$ that I would keep using for the rest. I fixed $lambda = 1$ and reconstructed the Shepp-Logan phantom for different values of $epsilon$, measuring the final error and the visual quality of the reconstruction.

#figure(
  grid(
    columns: (3fr, 1fr),
    align: horizon,
    inset: 10pt,
    image(
      "../figures/eps_error.svg",
      width: 75%
    ),
    table(
      columns: (auto, auto),
      inset: 10pt,
      align: center,
      table.header(
        [*$epsilon$*], [*Final Error*],
      ),
      [$0.001$], [$0.0402$],
      [$0.01$], [$0.0405$],
      [$0.1$], [$0.0528$],
      [$1.0$], [$0.0960$]
    ),
    grid.cell(
      colspan: 2,
      image(
        "../figures/eps_recos.svg",
      ),
    )

  ),
  
  caption: [Plot of the error vs. iterations, table of the final error, and reconstructions visualized for different values of $epsilon$ with $lambda = 1$]
)<fig:epsilon_choice>

On @fig:epsilon_choice, we can see that the final error is the lowest for $epsilon = 0.001$ and $epsilon = 0.01$. However, the reconstruction for $epsilon = 0.001$ looks smoother and less noisy than the one for $epsilon = 0.01$. Therefore, I chose to use $epsilon = 0.001$ for the rest of the project. One could have also argued that $epsilon = 0.01$ is a better choice because it converges much faster than $epsilon = 0.001$ but I wanted to prioritize the quality of the reconstruction over the speed of convergence.

Once $epsilon$ is fixed, I used Morozov's discrepancy principle to choose a good value for $lambda$. Since the noise level is known (because it was added artificially), I simply searched for the first greatest value of $lambda$ such that the residual is less than the noise level up to a certain tolerance $tau > 1$ which I arbitarily chose to be equal to $1.1$.

#figure(
 grid(
    rows: (auto, auto, auto),
    inset: 10pt,
    image(
      "../figures/morozov_2d.svg",
    ),
    table(
      rows: (auto, auto, auto),
      columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
      inset: 10pt,
      align: center,

      [*$lambda$*], [$10$], [4.6416], [$2.1544$], [#text(red)[$1$]], [$0.4642$], [$0.2154$], [$0.1$], [$0.0464$], [$0.0215$], [$0.01$],
      [*$tau$*], [$1.958$], [$1.393$], [$1.119$], [#text(red)[$1.015$]], [$0.955$], [$0.862$], [$0.781$], [$0.737$], [$0.717$], [$0.709$],
      [*Final Error*], [$0.1690$], [$0.1144$], [$0.0695$], [#text(red)[$0.0402$]], [$0.0277$], [$0.0369$], [$0.06$], [$0.0861$], [$0.1120$], [$0.1318$]

    ),

    image(
      "../figures/morozov_2d_reco.svg",
      width: 50%
    )

  ),
  
  caption: [Plots of the residual against $lambda$ in logscale, of the relative error against $lambda$ in logscale, and of the reconstructions visualized for different values of $lambda$ with $epsilon = 0.001$] 
)

Here, $lambda = 1$ was the first value such that the residual is less than the noise level up to a tolerance of $1.1$. Therefore, I chose to use $lambda = 1$ for the rest of the project. Notice that the final error is actually lowest for $lambda = 0.4642$ thus the method over-regularizes the solution. With $tau > 1$, we ensure that we stop on the safe side.

== FBP: weights and cutoff frequency

We were asked in the assignment to use a short-scan cone-beam geometry. This geometry covers $180 degree + 2 gamma_m$ where $gamma_m$ is the opening angle of the cone beam. This implies that some rays measure the same line integral of the object. However, the FBP method does not take this into account and thus we need to weight the rays to avoid over-counting. The weights are given by the Parker weights @parker-weights which are implemented by ODL. Their weights, however, are in flipped order. I therefore had to reimplement them to get the correct weights. 

#figure(
  image(
    "../figures/parker_check.svg",
    width: 75%
  ),
  caption: [Reconstruction using FBP without weighting, with ODL's Parker weights, with my implementation of the weights, and with a full $360 degree$ scan]
)<fig:parker_check>

On @fig:parker_check, we can see that the reconstruction without the weights is smoothed out (in the white region). Using ODL's weights, it actually worsens the reconstruction. On the other hand, using my implementation of the weights (which is just the ODL's parker weights but flipped), the reconstruction is much better. Finally, using a full $360 degree$ scan (without weighting), we can see that the reconstruction is as good as the short-scan with the correct weights.

Now that the we have the correct way of applying FBP, we can combine it with a low-pass filter to remove high-frequencies. The assignment asks to use Ram-Lak and Hann filters and compare them. To find a good cutoff frequency, I tried different values and measured the final error and the visual quality of the reconstruction which we can see on @fig:fbp_2d. We can see that the Hann filter gives significantly better results than the Ram-Lak filter. Even though the final error is lower for the 1.0 cutoff frequency, the reconstruction with 0.5 is a bit smoother.
Either way, the reconstructions are not as good as the ones obtained with the regularization method.

#figure(
  grid(
    rows: (auto, auto),
    image(
      "../figures/fbp_2d_Hann.svg",
      width: 75%
    ),
    image(
      "../figures/fbp_2d_Ram-Lak.svg",
      width: 75%
    ),
  ),
  caption: [Plots of the reconstructions visualized for different values of the cutoff frequency with Hann and Ram-Lak filters]
)<fig:fbp_2d>


#wrap-content(
  align: top + right,
  [#figure(
    image(
      "../figures/comparison_2d.svg",
      width: 100%
    ),
    caption: [Comparison of the reconstructions with the regularization method and with FBP]
  )<fig:comparison_2d>],



  [On @fig:comparison_2d, we can see that the regularization method gives a much better reconstruction than FBP. The regularization method recovers a much smoother image with the main ellipses clearly visible.]
)




= Reconstruction in 3D

Now, we can move on to the reconstruction of the 3D CT scan of a human head. The parameters $epsilon = 0.001$ was kept from the 2D case, same goes for the step size $t = 1.9/L$. The parameter $lambda$ was chosen using Morozov's discrepancy principle as before and thus need to be re-tuned for this case.

== Choice of $lambda$

As before, I used Morozov's discrepancy principle to choose a good value for $lambda$ but unlike the 2D case where a test image and artificial noise were used, here we only have access to the noisy sinogram. The assignment stated that the noise's relative norm is about $1%$. However, no matter the value of $lambda$, the residual stopped decreasing and was floored at a value of $approx 1.34 delta$. I thus ran the same procedure as before to obtain that floor value, considered the effective noise level to be $delta_"eff" = 1.34 delta$, and finally found the first value of $lambda$ such that the residual is less than $delta_"eff"$ up to a tolerance of $tau =1.1$.

On @fig:morozov_3d, we can see on the plot the initial noise estimate of $1%$ times the tolerance which the residual never reaches. On the other hand, the effective noise level of $1.34%$ times the tolerance is first reached for $lambda = 1.0$.

#figure(
  grid(
    rows: (auto, auto),
    inset: 10pt,
    table(
      rows: (auto, auto),
      columns: (auto, auto, auto, auto, auto, auto, auto, auto),
      inset: 10pt,
      align: center,

      [*$lambda$*], [$2.0$], [$1.5$], [#text(red)[$1.0$]], [$0.3$], [$0.1$], [$0.01$], [$0.001$],
      [*Residual / delta*], [$1.532$], [$1.488$], [#text(red)[$1.445$]], [$1.370$], [$1.347$], [$1.341$], [$1.340$],
    ),
    image(
      "../figures/morozov_3d.svg",
      width: 90%
    )
  ),
  caption: [Table of the residual against $lambda$ and plot of the residual against $lambda$ in logscale]
)<fig:morozov_3d>


== Comparison with FBP

We can now compare the reconstruction obtained with the regularization method with the one obtained with FBP with both filters previously used. To highlight the worse filter type and cutoff frequency, I chose to use the Ram-Lak filter with a cutoff frequency of $1$. I chose the Hann filter with a cutoff frequency of $0.25$ for comparison too because it gave the best results for the FBP method. The results can be seen on @fig:comparison_3d.

#figure(
  grid(
    rows: (auto, auto, auto),
    image(
      "../figures/final_fbp_ramlak.svg",
      width: 100%
    ),
    image(
      "../figures/final_fbp_hann.svg",
      width: 100%
    ),
    image(
      "../figures/final_huber.svg",
      width: 100%
    ),
  ),

  caption: [Comparison of the reconstructions with the regularization method and with FBP]
)<fig:comparison_3d>

Even though the FBP method with the Hann filter gives a better reconstruction than with the Ram-Lak filter and a very good reconstruction overall, the regularization method really smooths the uniform regions but still preserves the edges, giving this "cartoon"-like effect the assignment is telling about.


= Conclusion

In this project, we reconstructed a 3D CT scan of a human head from a noisy sinogram by minimizing a least-squares term regularized with the Huber function using gradient descent. The 2D Shepp-Logan phantom was first used to tune the method, which gave $epsilon = 0.001$ and a step size $t = 1.9/L$. In 3D, the residual never reached the given noise level, so I used an effective noise level $delta_"eff" = 1.34 delta$ in Morozov's discrepancy principle instead, which gave $lambda = 1$.

In both 2D and 3D, the regularization method gives better reconstructions than FBP. The downside is the computation time: FBP only takes a few seconds while the regularization method needs tens of minutes on the 3D volume. FBP can thus still be useful as a quick preview or as a starting point for the iterative method.

The main message of this project is that simply inverting the forward operator is not enough when the data is noisy. Adding what we know about the image we want to recover (here, mostly uniform regions separated by sharp edges) through a regularizer gives much better results.



#pagebreak()
#bibliography("refs.yml", full: true, title: "References")