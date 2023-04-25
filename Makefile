.PHONY: 

all :
	echo "Use make_ssl2.sh."

aprioricost :
	cd bblas_src_aprioricost ; \
	make -j8

benchmark : 
	cd benchmark ; \
	make -j8 

clean:
	-cd lib ; \
	rm -f *.o *.a *.so
	-cd bblas_src_aprioricost ; \
	rm -f *.o *.a *.so
	-cd benchmark ; \
	make clean

distclean:
	make clean
	-rm -rf bblas_src_aprioricost
