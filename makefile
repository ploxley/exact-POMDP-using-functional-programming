import: exact-pomdp-import.sml LP-solve-import.c
	mlton -default-ann 'allowFFI true' -export-header export.h  -link-opt -lglpk exact-pomdp-import.sml LP-solve-import.c
	
clean: 
	rm exact-pomdp-import output.dat
