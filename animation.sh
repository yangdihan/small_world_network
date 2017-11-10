ffmpeg -f image2 -r 5 -i *.png -vcodec mpeg4 -y movie.mp4
# rm *.png
python3 plot.py