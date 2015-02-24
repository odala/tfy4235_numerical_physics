PROGRAM dayNumber

	IMPLICIT NONE	!// Prevent the use of implicit declaration
	
	INTEGER					:: counter
	INTEGER, DIMENSION(12) 	:: months
	INTEGER					:: day, month
	INTEGER					:: daynr

	daynr = 0
	day = 16
	month = 9
	months(1) = 31
	months(2) = 28
	months(3) = 31
	months(4) = 30
	months(5) = 31
	months(6) = 30
	months(7) = 31
	months(8) = 31
	months(9) = 30
	months(10) = 31
	months(11) = 30
	months(12) = 31
	
	DO counter = 1, month - 1
		daynr = daynr + months(counter)
	END DO
	
	daynr = daynr + day
	
	PRINT *, daynr

END PROGRAM dayNumber
