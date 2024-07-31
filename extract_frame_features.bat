cif %1 -features all -header 1 > output\features_%2_frames_1.csv
for %%f in (10 20 40 80 160 320 640 1280) do cif %1 -features all -header 1 -frames %%f > output\features_%2_frames_%%f.csv