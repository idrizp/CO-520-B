
from PIL import Image
import glob
import os

out_dir = ''

# This converts all .BMP files in all directories to .png files. This is useful for hardcopy images from the oscilloscope.

for img in glob.glob('**/*.BMP', recursive=True):
    # print(img)
    path = os.path.join(os.getcwd(), str(img[:-4]) + '.png')
    print(path)
    Image.open(os.path.join(os.getcwd(), img)).save(path)

