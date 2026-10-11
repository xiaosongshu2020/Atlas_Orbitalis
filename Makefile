#===============================================================================
# Makefile for Atlas Project
# Compiler: gfortran
#===============================================================================

FC = gfortran
AR = ar
ARFLAGS = rcs
RM = rm -f
FFLAGS = -O2 -std=gnu -fall-intrinsics -fPIC -fopenmp

# 目录设置
SRC_DIR = src
TEST_DIR = test
APP_DIR = app
BUILD_DIR = build
BIN_DIR = bin

#===============================================================================
# 库模块编译（按依赖顺序）
#===============================================================================

LEVEL0 = module_const.f90 module_toolKits.f90 module_ErrorHandler.f90 \
         module_progressBar.f90 module_rotator.f90 module_randum.f90 \
         module_String.f90 module_FileDef.f90 module_elements.f90 \
         module_MathsMatrix.f90 module_sp3.f90 \
         module_EarthRadiationPressure.f90 module_shadow.f90

LEVEL1 = module_CRTBP.f90 module_lambert.f90 module_time.f90 module_EOPParam.f90

LEVEL2 = module_matrix.f90 module_planets.f90 module_halo.f90 \
         module_LissajousOrbit.f90 module_timeSeries.f90 module_clock.f90

LEVEL3 = module_space.f90 module_SLE.f90 module_attitude.f90 module_cartesian.f90

LEVEL4 = module_Gravity.f90 module_atmosphere.f90 \
         module_perturbationPN.f90 \
         module_perturbationSRP.f90 module_fastPrediction.f90 \
         module_perturbationThirdBody.f90 \
         module_thruster.f90 module_orbitPredictionNumerical.f90

LEVEL5 = module_station.f90 module_satellite.f90 module_access.f90 \
         module_residual.f90 module_satelliteScenario.f90 \
         module_IODor.f90 module_ODor.f90

LEVEL6 = module_Atlas.f90

FOR_FILES = PLEPH.for nrlmsise00_driver.for nrlmsise00_sub.for

ALL_SRC = $(LEVEL0) $(LEVEL1) $(LEVEL2) $(LEVEL3) $(LEVEL4) $(LEVEL5) $(LEVEL6) $(FOR_FILES)
OBJ_FILES = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(filter %.f90,$(ALL_SRC))) \
            $(patsubst %.for,$(BUILD_DIR)/%.o,$(filter %.for,$(ALL_SRC)))

ATLAS_LIB = $(BUILD_DIR)/libatlas.a

#===============================================================================
# SOFA 库（Fortran 77 源码在 sofa/ 下，249 个文件，每文件一个例程）
#
# 由源码编译成 build/libsofa.a，供所有可执行文件链接。它与 Atlas 相互独立、
# 不依赖任何 Atlas 模块，因此单独一条规则，先于程序构建完成。
# t_sofa_f.for 是 SOFA 自带的测试主程序，不编入库。
#===============================================================================

SOFA_DIR  = sofa
SOFA_SRCS = $(filter-out $(SOFA_DIR)/t_sofa_f.for,$(wildcard $(SOFA_DIR)/*.for))
SOFA_LIB  = $(BUILD_DIR)/libsofa.a
SOFA_OBJS = $(patsubst $(SOFA_DIR)/%.for,$(BUILD_DIR)/sofa/%.o,$(SOFA_SRCS))
SOFA_FLAGS = $(SOFA_LIB)

$(BUILD_DIR)/sofa:
	@mkdir -p $@

$(BUILD_DIR)/sofa/%.o: $(SOFA_DIR)/%.for | $(BUILD_DIR)/sofa
	$(FC) $(FFLAGS) -c $< -o $@

$(SOFA_LIB): $(SOFA_OBJS)
	$(AR) $(ARFLAGS) $@ $^

#===============================================================================
# 自动发现 app/ 和 test/ 下的程序
#===============================================================================

APP_SRCS = $(wildcard $(APP_DIR)/*.f90)
TEST_SRCS = $(wildcard $(TEST_DIR)/*.f90)

# 可执行文件后缀：MinGW 下 gfortran 会自动给 -o 指定的名字补上 .exe。
# 目标名必须与真实产物一致，否则规则声称产出 bin/xxx 而磁盘上只有 bin/xxx.exe，
# 该目标永远"不存在"，make 每次都会把它重新编译、重新链接一遍。
EXE = $(if $(findstring mingw,$(shell $(FC) -dumpmachine 2>/dev/null)),.exe,)

# 输出到 bin/ 目录
APP_BINS = $(patsubst $(APP_DIR)/%.f90,$(BIN_DIR)/%$(EXE),$(APP_SRCS))
TEST_BINS = $(patsubst $(TEST_DIR)/%.f90,$(BIN_DIR)/%$(EXE),$(TEST_SRCS))

#===============================================================================
# 目标
#===============================================================================

.PHONY: all clean

# 默认目标显式指定，避免被前面定义的规则（如 $(BUILD_DIR)/sofa）抢走
.DEFAULT_GOAL := all

all: $(BUILD_DIR) $(BIN_DIR) $(SOFA_LIB) $(ATLAS_LIB) $(APP_BINS) $(TEST_BINS)
	@echo "Build complete!"

$(BUILD_DIR) $(BIN_DIR):
	@mkdir -p $@

$(ATLAS_LIB): $(OBJ_FILES)
	$(AR) $(ARFLAGS) $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.for | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR)

$(BIN_DIR)/%$(EXE): $(APP_DIR)/%.f90 $(ATLAS_LIB) $(SOFA_LIB) | $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $(BUILD_DIR)/$*.o
	$(FC) $(FFLAGS) $(BUILD_DIR)/$*.o $(ATLAS_LIB) $(SOFA_FLAGS) -o $@

$(BIN_DIR)/%$(EXE): $(TEST_DIR)/%.f90 $(ATLAS_LIB) $(SOFA_LIB) | $(BIN_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) -I$(BUILD_DIR) -c $< -o $(BUILD_DIR)/$*.o
	$(FC) $(FFLAGS) $(BUILD_DIR)/$*.o $(ATLAS_LIB) $(SOFA_FLAGS) -o $@

ifneq ($(EXE),)
# 便利规则：make bin/test_CRTBP 等价于 make bin/test_CRTBP.exe（必须放在上面两条
# 带后缀的规则之后，make 按定义顺序选取首个前置可满足的模式规则）
$(BIN_DIR)/%: $(BIN_DIR)/%$(EXE)
	@:
endif

clean:
	$(RM) -rf $(BUILD_DIR) $(BIN_DIR)

help:
	@echo "Usage: make        - Build all"
	@echo "       make clean  - Remove build directories"
