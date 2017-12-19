################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Sources/Hermitian3x3.cpp \
../Sources/Matrix3x3.cpp \
../Sources/main.cpp 

OBJS += \
./Sources/Hermitian3x3.o \
./Sources/Matrix3x3.o \
./Sources/main.o 

CPP_DEPS += \
./Sources/Hermitian3x3.d \
./Sources/Matrix3x3.d \
./Sources/main.d 


# Each subdirectory must supply rules for building sources it contributes
Sources/%.o: ../Sources/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0  -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


