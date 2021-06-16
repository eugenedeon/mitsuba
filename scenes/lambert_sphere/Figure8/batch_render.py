import sys
import os
import math

intensity = 250

kdred = 0.5
kdgreen = 0.25
kdblue = 0.12

count = 0
for rot in [0,30,60,90,105,120,135,150,160]:
  os.system( 'mitsuba -o lambert_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' lambert.xml' )
  os.system( 'mitsuba -o diffusionT_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' diffusion_transport.xml' )
  os.system( 'mitsuba -o lambert_sphere_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' lambert_sphere.xml' )
  os.system( 'mitsuba -o lambert_sphere_fast_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' lambert_sphere_fast.xml' )
  os.system( 'mitsuba -o lambert_sphere_hapke_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' lambert_sphere_hapke.xml' )
  os.system( 'mitsuba -o ON_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' ON.xml' )
  os.system( 'mitsuba -o chandra_' + str(count) + '.jpg -Drot=' + str(rot) + ' -Dkdred=' + str(kdred) + ' -Dkdgreen=' + str(kdgreen) + ' -Dkdblue=' + str(kdblue) + ' -Dintensity=' + str(intensity) + ' chandra.xml' )

  count += 1
  