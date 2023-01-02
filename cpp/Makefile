diffusion: diffusion.cpp
	g++ diffusion.cpp -L/usr/local/lib/ -lsundials_nvecserial -lsundials_ida \
    -o diffusion

jac: jac.cpp
	g++ jac.cpp -L/usr/local/lib/ -lsundials_nvecserial -lsundials_ida \
    -o jac

test: test_band.cpp
	g++ test_band.cpp -L/usr/local/lib/ -lsundials_nvecserial -lsundials_sunmatrixband \
    -lsundials_sunlinsolband -o test
