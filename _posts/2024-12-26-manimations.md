---
layout: post
math: true
title:  "Manim, and variations of a geometry problem"
date:   '2024-12-26'
categories: python manim geometry
---

## IMC 2024 problem 5

Every year the University College London hosts the [International Mathematics
Competitions for College Students (IMC)](https://www.imc-math.org.uk), which is
often hosted in Blagoevgrad, Bulgaria. I participated back in 2013 and, while I
am not at all proud of my result, I still like to occasionally check the recent
exams to see if I can solve any problem (I almost always cannot).

The competition is hosted in two days, each one of which students have 5 hours
to solve 5 problems, each worth 10 points. Your total score is the sum of points
you get in each problem, out of 100. The first problem is usually the easiest,
and the fifth problem the hardest. This is obviously relative but statistically
holds out in almost every year. In 2024, [problem 5 of day
1](https://www.imc-math.org.uk/?year=2024&section=problems&item=day1), the
"hardest" of the first day, reads:


<span style="color:#DDDDDD">
Let $$n > d$$ be positive integers. Choose $$n$$ independent, uniformly
distributed random points in $$x_1, \dots, x_n$$ in the unit ball $$B \subset
\mathbb{R}^d$$ centered at the origin. For a point $$p \in B$$, denote by
$$f(p)$$ be the probability that the convex hull of $$x_1, \dots, x_n$$ contains
$$p$$. Prove that, if $$p,q \in B$$ and the distance of $$p$$ from the origin is
smaller than the distance of $$q$$ from the origin, then $$f(p) \geq f(q)$$.
</span>


Putting it simply, it means that, in any arbitrary dimension $$d$$, the closer a
point is to the origin, the more likely it is to be contained in a convex hull
of random points sampled from the unit ball (the condition $$n > d$$ guarantees
that the convex hull is also $$d$$-dimensional).

Even after reading the
[solution](https://www.imc-math.org.uk/?year=2024&section=problems&item=prob5s),
I don't think I can intuitively understand why. Usually when we think about
these problems it helps to set the smallest possible values for the variables,
so let's think about the case where $$d=2$$, $$n=3$$. Obviously, if $$p$$ is in the
border of the circle, the probability is $$0$$, and a corollary of the problem
statement is that the origin ($$p=(0,0)$$) would be the most likely point to be
contained in the convex hull.

## Simplifying the problem

Let's take a quick detour and try to find $$f(p)$$ when $$p=0$$ in this case,
in other words, the probability that a triangle from three points in a unit
circle contains the origin. These three points can be defined by angles
$$\theta_1$$, $$\theta_2$$ and $$\theta_3$$ such that coordinates $$(x_i,
y_i)$$ are $$(\mathrm{cos}\theta_i, \mathrm{sin}\theta_i)$$, for $$i \in
\{1,2,3\}$$. Without loss of generality we can assume $$\theta_3 = 0$$ and
$$\theta_2 > \theta_1$$, so that we are picking two angles that split the
circle in three arcs. It should be easy to see that the origin is inside the
triangle if it is acute, which means that the angles for the three arcs the
triangle splits must be no more than $$\pi$$, as shown below, where the left
triangle does not contain the origin because $$\theta_3 > \pi$$.

![The left triangle does not contain the origin because $$\theta_3 > \pi$. The right triangle contains the origin because the angle in all arcs is no more than $$\pi$$.](https://i.ibb.co/b7C7LBP/example.png)



Mathematically, if $$\theta_i$$ are uniformly sampled angles that add to $$2
\pi$$, it means that $$\theta_1 - \theta_3 < \pi$$, $$\theta_2 - \theta_1 <
\pi$$ and $$2\pi - \theta_2 < \pi$$. Since $$\theta_3 = 0$$ this is reduced to
$$\theta_1 < \pi$$ and $$\pi < \theta_2 < \pi + \theta_1$$, with $$\theta_2
\sim \mathrm{Uniform}(0, 2\pi)$$, so integrating over all values of
$$\theta_1$$, our answer is

$$
2 \times \frac{1}{2\pi}\int_{0}^{\pi} \frac{\theta}{2\pi} d\theta = 2 \times \frac{\pi^2}{2} \times \frac{1}{4\pi^2} = \frac{1}{4},
$$

where the first multiplication by $$2$$ accounts for the symmetry condition
where $$\theta_2 < \theta_1$$, which occurs with equal probability as $$\theta_1
< \theta_2$$. So according to the corollary of the problem, $$f(p)$$ is a
decreasing function of $$\|p\|$$ from $$f(p)=0$$ if $$\|p\|=1$$ and $$f(p)=1/4$$
when $$\|p\|=0$$. As should be expected, this argument expands to higher
dimensions, where we can generalize that in dimension $$d$$, the probability
that a simplex with $$d+1$$ vertices on the surface of a unit ball in
$$\mathbb{R}^d$$ is $$2^{-d}$$. For example, for two dimensions, $$d = d$$ and
$$2^{-2}=1/4$$ as discussed. You can read more about the proof
[here](https://mathworld.wolfram.com/SphereTetrahedronPicking.html).

As much as I'd like to find the formula for $$f(p)$$ for arbitrary values of
$$n$$ and $$d$$, I'm afraid this is as far as my capacity to think about this
problem goes.

Now onto the real purpose of this post.

## Manim

Manim is a mathematical animation package for Python which facilitates the
creation of educational videos for mathematics. It is used and co-developed by
the amazing youtube creator of mathematical videos
[3blue1brown](https://www.youtube.com/c/3blue1brown), and used by several other
YouTube channels. I always wanted to learn how to use it, and I will use the
above problem to create a basic animation with the [manim python
package](https://docs.manim.community/en/stable/index.html). Here is our goal:

**Animate 50 uniformly random triangles in an unit circle, color them green if
they contain the origin, and red if not**.


I installed manim using [conda](https://anaconda.org/conda-forge/manim), which
includes both the binary to compile the code and the python packages.

After studying the quickstart basics, I arrived at the following code.

```python
# practice.py
# author: Guilherme Sena

from manim import *
import numpy as np
import random

# A radius that covers the entire video
RADIUS=6

# convert to polar coordinates
def polar_point(theta, radius=RADIUS):
    return [radius*np.cos(theta), radius*np.sin(theta), 0]

# makes three triangle points with polar coordinates
def polar_triangle(t1, t2, t3, radius=RADIUS):
    p1 = polar_point(t1, radius)
    p2 = polar_point(t2, radius)
    p3 = polar_point(t3, radius)
    
    return [p1, p2, p3]

# checks if the arcs are all below PI, meaning that the triangle
# contains the origin
def has_origin(p):
    a1 = p[1] - p[0]
    a2 = p[2] - p[1]
    a3 = 2*PI - a1 - a2
    return a1 <= PI and a2 <= PI and a3 <= PI

# creates IID random angles in increasing order
def random_angles():
    v = [
        random.uniform(0, 2*PI),
        random.uniform(0, 2*PI),
        random.uniform(0, 2*PI)
    ]
    v.sort()
    return v

# creates a random manim triangle with the above functions
def random_triangle():
    rand_angles = random_angles()
    
    p = polar_triangle(rand_angles[0], rand_angles[1], rand_angles[2])
    triangle = Polygon(*p)
    
    c = GREEN if has_origin(rand_angles) else RED
    triangle.set_fill(color=c, opacity=0.5)
    triangle.set_stroke(color=c, width=3)   

    return triangle    

# manim animations need to override the "Scene" class, this is the
# only class in this python file, so we don't need to pass "manim Practice"
# when running
class Practice(Scene):
    def construct(self):
        circle = Circle(radius=RADIUS, color=BLUE, fill_opacity=0.7)
        polygon = random_polygon()
        circle.move_to(ORIGIN)
        
        self.add(circle)
        self.add(triangle)
        self.add(Dot(circle.get_center()).set_color(YELLOW).scale(1.5))

        for i in range(50):
            t = random_triangle()
            self.play(Transform(triangle, t))
            self.wait(0.1)
```

This is how you compile to create a `512 x 512` GIF with the animation.

```bash
$ manim -i -r 512,512 practice.py
```

And this is the final result.

![Random triangles in an unit circle](https://i.ibb.co/6s0vTn3/Practice-Manim-CE-v0-18-1.gif)

As an addition, using the functions above I made a "computational" proof that
the probability of the origin being inside the triangle is, indeed, $$1/4$$:

```python
s = 0
lim = 1000000
for i in range(lim):
    if has_origin(random_angles()):
        s  += 1
        
print("%s %s %s" % (s, lim, s/lim))
# 249858 1000000 0.249858
```

My next step is to animate in such way that the triangle rotates its points
around the corner of the circle instead of a linear interpolation.


