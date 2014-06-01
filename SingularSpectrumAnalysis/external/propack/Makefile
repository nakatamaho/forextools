PROPACKVER=2.1
BUILDDIR=PROPACK

EXTRACT_STAMP=.propack_extract_done
PATCH_STAMP=.propack_patch_done
CONFIGURE_STAMP=.propack_configure_done
BUILD_STAMP=.propack_build_done
INSTALL_STAMP=.propack_install_done
CHECK_STAMP=.propack_check_done

all: extract config build check install

work/$(EXTRACT_STAMP):
	rm -rf work
	mkdir -p work
	cd work; tar xvfz ../../archives/PROPACK77_v${PROPACKVER}.tar.gz
	touch work/$(EXTRACT_STAMP)

work/$(PATCH_STAMP): work/$(EXTRACT_STAMP)
	cd work/$(BUILDDIR);  patch -p0 < ../../../archives/patch-propack
	touch work/$(PATCH_STAMP)

work/$(CONFIGURE_STAMP): work/$(PATCH_STAMP)
	cd work/$(BUILDDIR); cp Make/make.linux_gcc_ia32 make.inc
	touch work/$(CONFIGURE_STAMP)

work/$(BUILD_STAMP): work/$(CONFIGURE_STAMP)
	cd work/$(BUILDDIR); make
	touch work/$(BUILD_STAMP)

work/$(CHECK_STAMP): work/$(BUILD_STAMP)
	cd work/$(BUILDDIR) ; make test ; make verify
	touch work/$(BUILD_STAMP)

work/$(INSTALL_STAMP): work/$(CHECK_STAMP)
	cd work/$(BUILDDIR); cp double/libdpropack_MINGW32.a ../../../../
	touch work/$(INSTALL_STAMP)

extract: work/$(EXTRACT_STAMP)

patch: work/$(PATCH_STAMP)

config: work/$(CONFIGURE_STAMP)

build: work/$(BUILD_STAMP)

check: work/$(CHECK_STAMP)

install: work/$(INSTALL_STAMP)

clean:
	rm -rf work

distclean:
	rm -rf work

