
set(CMAKE_BUILD_TYPE Debug)

# The test files I care only (create a list)
set(testsources
	test_init.c
	test_set_strain.c
	test_bc_apply_on_u.c
	test_bc_apply_on_u_omp.c
	test_assembly_res.c
	test_assembly_res_omp.c
	test_assembly_jac.c
	test_newton_raphson.c
	test_newton_raphson_omp.c
	test_calc_stress.c
	test_calc_strain.c
	test_solver_linear.c
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
add_test(NAME test_assembly_res COMMAND test_assembly_res)
add_test(NAME test_assembly_jac COMMAND test_assembly_jac)
add_test(NAME test_newton_raphson COMMAND test_newton_raphson)
add_test(NAME test_calc_stress COMMAND test_calc_stress)
add_test(NAME test_calc_strain COMMAND test_calc_strain)
add_test(NAME test_solver_linear COMMAND test_solver_linear 10)
