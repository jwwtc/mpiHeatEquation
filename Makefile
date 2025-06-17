CC?=mpicc
CFLAGS=-O2 -std=c99
SRC=src/main.c src/functions.c
TARGET=main.exe
SERIAL_SRC=extra/ser_heat_equation.c
SERIAL_TARGET=ser_heat.exe

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $(SRC) -lm

clean:
	rm -f $(TARGET) $(SERIAL_TARGET)

serial: $(SERIAL_TARGET)

$(SERIAL_TARGET): $(SERIAL_SRC)
	gcc $(CFLAGS) -o $@ $(SERIAL_SRC) -lm

.PHONY: all clean serial
