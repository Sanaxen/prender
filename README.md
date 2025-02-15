# prender
<center><img src=./images/image00.png width=80%></center>

prender is an impractical renderer that I created while experimenting with various ways to learn how to implement path tracing.

Very much learned from the following sites
[https://kagamin.net/hole/edupt/](https://kagamin.net/hole/edupt/)
This site can handle diffuse reflections, specular surfaces, and ideal glass.
[https://github.com/githole/edupt[(https://github.com/githole/edupt)]


## prender
Prender adds subsurface scattering, translucency, and involved media to these.
The light source can also handle collimated light sources and spotlights.
In addition, experimental implementations can handle full-spectrum rendering, gravity rendering (black holes), and wormholes.

Gravity rendering (black holes) and wormholes are influenced by the movie Interstellar.
The following papers are used as references for implementation.  
[Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar](https://arxiv.org/abs/1502.03808)  

[Visualizing Interstellar's Wormhole](https://arxiv.org/abs/1502.03809)

**Gravitational Renderer**  `The VFX used in Interstellar was a renderer (Double Negative Gravitational Renderer) developed by 30 members of Double Negative's VFX team under Kip's guidance.
 The technology developed for this film and the new findings that resulted in two papers also attracted a lot of attention.
 Gravitational Lensing by Spinning Black Holes in Astrophysics, and in the Movie Interstellar and Visualizing_Interstellar's_Wormhole.
 It describes how a spinning black hole or wormhole looks when approached (and what kind of image the camera records).
 Director Christopher Nolan discussed with Kip how to represent the vicinity of a black hole in VFX,
  i.e., how to depict space-time in a strong gravitational field, and eventually decided to develop the above renderer to be faithful to physics.

As for wormholes, I was very inspired by the following  
[Space Time Travel](https://www.spacetimetravel.org/wurmlochflug)
[Movie: Flight through a wormhole](https://www.youtube.com/watch?v=SZDOKtT_QZE)


[Wormhole Simulation](https://www.youtube.com/embed/SZDOKtT_QZE)

---
# **[prender examples](https://github.com/Sanaxen/prender/blob/main/example.md)**

## Notes.
prender was developed around 2014.
At that time, we were working on various projects and trials, so there may be some incomplete parts or code based on incorrect ideas.

## There are some things we really need to be careful about.  
It has absolutely no practical value in the following respects.
1. execution speed is very slow.
There are experimental implementations of MPI and other mechanisms, but they are basically just parallelization with multi-threading. A state-of-the-art implementation would probably use GPUs, but that has not been done at all.
2. There is a lot of experimental thinking in the implementation.
Some parts of the theory are not implemented as they are, and some implementations are based on understanding supplemented by imagination.