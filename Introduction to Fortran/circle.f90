PROGRAM circle

	IMPLICIT NONE	!// Prevent the use of implicit declaration
	
	REAL				:: radius
	CHARACTER(LEN=42)	:: radius_prompt
	
	radius_prompt = "Enter the radius of your circle / sphere: "
	WRITE(*, FMT = '(A)', ADVANCE="NO") radius_prompt
	READ(*,*) radius
	
	PRINT *, "The circumference of the circle is: ", circumference_of_circle(radius)
	PRINT *, "The area of the circle is: ", area_of_circle(radius)
	PRINT *, "The volume of the sphere is: ", volume_of_sphere(radius)
	
CONTAINS

FUNCTION circumference_of_circle(radius) RESULT(circumference)
	IMPLICIT NONE
	REAL, INTENT(IN)		:: radius
	REAL					:: circumference
	! Declare local constant Pi
	REAL, PARAMETER 		:: pi = 3.1415927
	circumference = 2 * pi * radius
END FUNCTION circumference_of_circle

FUNCTION area_of_circle(radius) RESULT(area)
	IMPLICIT NONE
	REAL, INTENT(IN)		:: radius
	REAL					:: area
	! Declare local constant Pi
	REAL, PARAMETER 		:: pi = 3.1415927
	area = pi * radius * radius
END FUNCTION area_of_circle

FUNCTION volume_of_sphere(radius) RESULT(volume)
	IMPLICIT NONE
	REAL, INTENT(IN)		:: radius
	REAL					:: volume
	! Declare local constant Pi
	REAL, PARAMETER 		:: pi = 3.1415927
	volume = 4 * pi * radius * radius * radius * 3
END FUNCTION volume_of_sphere

END PROGRAM circle
