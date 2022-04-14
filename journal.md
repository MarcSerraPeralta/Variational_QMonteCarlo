# Weekly progress journal

## Instructions

In this journal you will document your progress of the project, making use of the weekly milestones.

Every week you should 

1. write down **on the day of the lecture** a short plan (bullet list is sufficient) of how you want to 
   reach the weekly milestones. Think about how to distribute work in the group, 
   what pieces of code functionality need to be implemented.
2. write about your progress **until Monday, 23:59** before the next lecture with respect to the milestones.
   Substantiate your progress with links to code, pictures or test results. Reflect on the
   relation to your original plan.

We will give feedback on your progress on Tuesday before the following lecture. Consult the 
[grading scheme](https://computationalphysics.quantumtinkerer.tudelft.nl/proj1-moldyn-grading/) 
for details how the journal enters your grade.

Note that the file format of the journal is *markdown*. This is a flexible and easy method of 
converting text to HTML. 
Documentation of the syntax of markdown can be found 
[here](https://docs.gitlab.com/ee/user/markdown.html#gfm-extends-standard-markdown). 
You will find how to include [links](https://docs.gitlab.com/ee/user/markdown.html#links) and 
[images](https://docs.gitlab.com/ee/user/markdown.html#images) particularly.

## Week 1

### Bullet list:

1. Research on the process of VMC and the different algorithms (@mserraperalta, @abermejillo, @dbedialaunetar)
6. Write a sketch for the structure of the library (@mserraperalta)
2. Search for a trial wavefunction and compute the analytic $`E_{local}`$ for the harmonic oscillator (@mserraperalta)
3. Implement the Metropolis algorithm (@abermejillo)
4. Implement the Monte Carlo integration algorithm (@dbedialaunetar)
5. Validate the results for E(parameter) for the harmonic oscillator (@mserraperalta, @abermejillo, @dbedialaunetar)

### Progress:

First, we will briefly comment how the tasks from the bulletlist have been performed and by whom.

1. @dbedialaunetar, @abermejillo and @mserraperalta all did some research and discussed how to best implement the different algorithms. 
2. @mserraperalta prepared a [skeleton](bc84c58e966ce216faaf9866dcb1db6ba76a8eba) for the project.
2. @mserraperalta implemented [symbolic calculation](cdf8cd38988a0217d29d56d02a2e6b06d675db6b) of the local E and the trial wavefunction. First he included the case of the Harmonic Oscillator.
3. @abermejillo implemented the [Metropolis algorithm](78626b745c015c17be106f9ed53da549e6d832f7) and also a function [find optimal trial move](b2b089f2ad6283aa2bcb87b90826f0e14d3e6175) that finds a trial move with average acceptance ratio 0.5 $`\pm`$ tol.
4. @dbedialaunetar implemented the [Monte Carlo integration algorithm](46fac493d8c56504cbc0633f53fa38f2289668e0).
5. @abermejillo and @mserraperalta checked that the walkers [sampled a gaussian](1794ed020ce402570b8f041dc9f4cbb914da3cf0) correctly and @mserraperalta that the energy of the [Harmonic Oscillator](f76f056e789b16f8b5173b34f587eee51aeee0ca) is computed correctly.

**Results and comments**

After implementing the Metropolis algorithm we have to chech that functions are properly sampled. We perform this analysis for a very simple function, a gaussian. The parameter that needs to be fixed is the trial move standard deviation, because we made a choice to perform trial moves according to a gaussian distribution, which is clearly balanced. 

$`\sigma_{tm}=10`$             |  $`\sigma_{tm}=1`$ 
:-------------------------:|:-------------------------:
![alt text](results/W1_walker_sig10.png) |  ![alt text](results/W1_walker_sig1.png)

$`\sigma_{tm}=0.01`$             |  $`\sigma_{tm,opt}=1.987`$ 
:-------------------------:|:-------------------------:
![alt text](results/W1_walker_sig01.png) |  ![alt text](results/W1_walker_sigopt.png)

In these figures we can see how first we sample a gaussian doing 50000 steps for three different values of the trial move. According to the choice of units, a distance of the order of unity if the most appropriate, and this is confirmed by the sampling. When the trial move is 10 the sampling is a bit worse than when it is 1, and the sampling for 0.01 is completely inservible. Finally, we compute an optimal trial move $`\sigma_{tm}=1.987`$, for which the sampling is adapted very nicely to the gaussian. 



(due 18 April 2022, 23:59)


## Week 2
(due 25 April 2022, 23:59)


## Week 3
(due 2 May 2022, 23:59)

