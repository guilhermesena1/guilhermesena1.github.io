---
layout: post
math: true
title:  "Setting up a blog"
date:   '2019-12-11'
categories: hello
---

I am configuring this blog if I want to communicate, write or document things
using markdown and figures. This post is to test if I can write some code with
syntax highlight and some math.

Here's an example c++ code:
~~~ cpp
// hi.cpp
#include <iostream>
int
main(int argc, const char **argv) {
  // prints hi to stdout
  std::cout << "Your program has " << argc << " arguments" << std::endl;
  return 0;
}
~~~

This is how to compile it:
~~~ bash
$ g++ -O3 -Wall -o hi hi.cpp
~~~

And this is how you run it
~~~ bash
$ ./hi arg1 arg2 arg3
Your program has 4 arguments
~~~

This is a [path to a link][path-to-link], and below is an equation:

$$
x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
$$

[path-to-link]: https://sena.works
