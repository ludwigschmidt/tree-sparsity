CXX=g++-5
CXXFLAGS=-std=gnu++11 -DNDEBUG -Wall -Wextra -march=native -O3 -I external/eigen -I src/include
GTEST_DIR=~/code/googletest/googletest

CORE_DIR=src/core
PYTHON_DIR=src/python
TEST_DIR=src/test
TEST_BIN_DIR=test_bin


clean:
	rm -rf obj
	rm -rf $(TEST_BIN_DIR)



obj/gtest-all.o: $(GTEST_DIR)/src/gtest-all.cc
	  mkdir -p obj
		  $(CXX) $(CXXFLAGS) -I $(GTEST_DIR)/include -I $(GTEST_DIR) -c -o $@ $<

obj/gtest_main.o: $(GTEST_DIR)/src/gtest_main.cc
	  mkdir -p obj
		  $(CXX) $(CXXFLAGS) -I $(GTEST_DIR)/include -c -o $@ $<

TEST_DEPS=$(TEST_DIR)/test_helpers.h obj/gtest-all.o obj/gtest_main.o

$(TEST_BIN_DIR)/exact_tree_projection_test: $(TEST_DIR)/exact_tree_projection_test.cc $(CORE_DIR)/exact_tree_projection.h $(TEST_DEPS)
	mkdir -p $(TEST_BIN_DIR)
	$(CXX) $(CXXFLAGS) -I $(GTEST_DIR)/include -c -o obj/exact_tree_projection_test.o $(TEST_DIR)/exact_tree_projection_test.cc
	$(CXX) $(CXXFLAGS) -o $(TEST_BIN_DIR)/exact_tree_projection_test obj/gtest_main.o obj/gtest-all.o obj/exact_tree_projection_test.o -pthread
