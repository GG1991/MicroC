
set(CMAKE_BUILD_TYPE Debug)

# The test files I care only (create a list)
set(testsources
	test_init.c
	test_set_strain.c
	test_bc_apply_on_u.c
	)

# Iterate over the list above
foreach (testfile ${testsources})
	# Delete File extensions (test_i.cpp -> test_i)
	get_filename_component(testname ${testfile} NAME_WE)

	# Create the executable test_i and link.
	add_executable(${testname} ${testfile})
	target_link_libraries(${testname} microc)
endforeach ()

# As some tests will require commands and could be executed in loops we add the
# tests individually here.

# Add a test that requires arguments
add_test(NAME test_init COMMAND test_init)
add_test(NAME test_set_strain COMMAND test_set_strain)
add_test(NAME test_bc_apply_on_u COMMAND test_bc_apply_on_u)
