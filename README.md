plot.py:
	input: the output .txt with "...false..." in name
	output: a plot of force vs. displacement diagram, with max-force labeled

./msh folder:
	contains various .msh

random_long_link.py:
	input: origin .msh and number of random links to add
	output: a modified .msh

add_link_tracker.txt:
average_x_before_add_long 	average_L_before_add_long
add_1_x_1		add_1_y_1		add_1_x_2		add_1_y_2
add_2_x_1		add_2_y_1		add_2_x_2		add_2_y_2
add_3_x_1		add_3_y_1		add_3_x_2		add_3_y_2
...				...				...				...



add_long_link_info.txt:
add_1_force		add_1_s		add_2_force		add_2_s		add_2_force		add_2_s
... every iteration ...