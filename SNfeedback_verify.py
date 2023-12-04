import os
import cv2 as cv


def within_cube(posx, posy, posz, cubexyz):
    if posx > cubexyz[0] and posx < cubexyz[1] and posy > cubexyz[2] and posy < cubexyz[3] and posz > cubexyz[4] and posz < cubexyz[5]:
        return True
    return False

def calc_potential_cube(target_posx, target_posy, target_posz):
    cubex_min = target_posx - 15
    cubex_max = target_posx + 15
    cubey_min = target_posy - 15
    cubey_max = target_posy + 15
    cubez_min = target_posz - 15
    cubez_max = target_posz + 15
    return [cubex_min, cubex_max, cubey_min, cubey_max, cubez_min, cubez_max]


# find a blob in images
# record the center coordinate and timestamp
target_posx = 0
target_posy = 0
target_posz = 0
timestamp = 200

# retrieve the cooresponding mask
mask_root = ""
mask_filename = ""
mask_file_path = os.path.join(mask_root, mask_filename)
mask = cv.imread(mask_file_path)

# form a proposal cube
cubexyz = calc_potential_cube(target_posx, target_posy, target_posz)


# search within SNfeedback.dat for SN center point within the cube range
SNfeedback_file = ""
SN_potential = []

with open(SNfeedback_file, "r") as f:
    data = f.readlines()

# remove the header from the data
data.pop(0)
for line in data:
    # if the center lies within the cube, record it
    line = line.split('\t')
    posx, posy, posz = line[5], line[6], line[7]

    # TODO: confirm if we need to consider timestamp as well
    if(within_cube(posx, posy, posz, cubexyz)):
        # record this SN
        SN_potential.append(line)


print(SN_potential)



