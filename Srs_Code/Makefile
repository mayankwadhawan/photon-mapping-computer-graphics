TARGET = bin_glut

SRC = \
	main.cpp \

INCLUDE = \
	-I./ \

LIB = 	\


LIBDEF= 	\
	-L/usr/X11R6/lib \
	-lglut -lGLU -lGL

DEFINE = \
	-D_LINUX_ \

CC = g++ -g

CFLAGS = \
	-O1 \
	-msse2\
	-g \
	-Wall \
	-funroll-all-loops\
	-ftree-vectorize\


OBJ = $(patsubst %.cpp,%.o,$(filter %.cpp,$(SRC)))

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(CFLAGS) $(DEFINE) $(INCLUDE) -c $< -o $@

all: $(TARGET)

libraries:
	@ for d in $(dir $(LIB)); do \
		make -C $$d; \
	done

$(TARGET): $(OBJ) $(LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB) $(LIBDEF)

clean:
	rm -f $(TARGET) $(OBJ)

clobber: clean
	@ for d in $(dir $(LIB)); do \
		make -C $$d clean; \
	done



