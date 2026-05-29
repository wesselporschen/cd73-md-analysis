load 4h2i_medoid_cluster_4.pdb
load hull_frame_24281.pdb
color cyan, frame_24281 and chain A
color green, frame_24281 and chain B
color hotpink, hull_frame_24281
select residue_selection, (resi 366 and chain A) or (resi 456 and chain B) or (resi 543 and chain A) or (resi 543 and chain B) or (resi 484 and chain A) or (resi 484 and chain B) or (resi 547 and chain B) or (resi 473 and chain A) or (resi 473 and chain B) or (resi 545 and chain A) or (resi 545 and chain B) or (resi 471 and chain A) or (resi 455 and chain A) or (resi 542 and chain A) or (resi 453 and chain A) or (resi 453 and chain B) or (resi 454 and chain A) or (resi 454 and chain B) or (resi 544 and chain A) or (resi 544 and chain B) or (resi 366 and chain A) or (resi 531 and chain A)
show sticks, residue_selection
color atomic, residue_selection
set stick_color, hotpink, residue_selection and elem C
zoom hull_frame_24281