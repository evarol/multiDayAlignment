# multiDayAlignment

Alignment of multiple day microscopy images.

![Demo](https://github.com/evarol/multiDayAlignment/blob/master/demo_2.png)


Step 1: Open multidayAlignment.m in Matlab, press run and select multiple images to align (.tif, .png files). It's important to select files in order you'd like them to be registered e.g. day 1 then day 2 then day 3 (Don't just multi select them in one shot).

Step 2: Press run and follow on screen instructions to find pairwise landmarks that align. When done press right mouse click (If you have no landmarks still press right click).

Step 3: When done with aligning all pairs of days run composeRegistrations.m. Select multiple images to align in the same order as you selected them in step 1.
