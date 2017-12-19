################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Tests/Hermitian3x3_test.cpp \
../Tests/Matrix3x3_test.cpp 

OBJS += \
./Tests/Hermitian3x3_test.o \
./Tests/Matrix3x3_test.o 

CPP_DEPS += \
./Tests/Hermitian3x3_test.d \
./Tests/Matrix3x3_test.d 


# Each subdirectory must supply rules for building sources it contributes
Tests/%.o: ../Tests/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0  -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


