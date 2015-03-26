#! /usr/bin/perl
#
# Usage: makemake {<program name> {<F90 compiler or fc or f77 or cc or c>}}
#
# Generate a Makefile from the sources in the current directory.  The source
# files may be in either C, CC, FORTRAN 77, Fortran 90 or some combination of
# these languages.  If the F90 compiler specified is cray or parasoft, then
# the Makefile generated will conform to the conventions of these compilers.
# To run makemake, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Orriginal version by Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# Modified by Ingve Simonsen <ingves@phys.ntnu.no> August, 2000
# NTNU, Trondheim, Norway 
#
# Mar 2007 : Added command to suppress the modula2 default.....
#


# Checks if a Makefile already exists
if (-e "Makefile")
{
		print "A Makefile already exists!   Overwrite (y/n) ";
		$answer = <STDIN>;
		if ( $answer =~ /^[nN]/ )
		{
		print "Exiting --- Makefile left unchanged\n";
		exit 0;
		}   
}


# Lets make the Makefile !
open(MAKEFILE, "> Makefile");
#
print MAKEFILE "PROG =\t$ARGV[0]\n\n";
#
# Source listing
#
print MAKEFILE "SRCS =\t";
@srcs = <*.f90 *.f *.F *.c *.cc *.cpp>;
&PrintWords(8, 0, @srcs);
print MAKEFILE "\n\n";
#
# Object listing
#
print MAKEFILE "OBJS =\t";
@objs = @srcs;
foreach (@objs) { s/\.[^.]+$/.o/ };
&PrintWords(8, 0, @objs);
print MAKEFILE "\n\n";
#

#----  C -------
print MAKEFILE "C           = gcc\n";
print MAKEFILE "CFLAGS      = -O\n";
print MAKEFILE "CLIBS       = \n";
print MAKEFILE "CLDFLAGS    = -s -L\$(HOME)/C/Lib/\n\n";
#----  CPP -----
print MAKEFILE "CPP         = g++\n";
print MAKEFILE "CPPFLAGS    = -O -I. -I\$(HOME)/CC/include/\n";
print MAKEFILE "CPPLIBS     = \n";
print MAKEFILE "CPPLDFLAGS  = -s -L\$(HOME)/CC/Lib/\n\n";
#----  F77 -----
print MAKEFILE "F77         = pgf77\n";
print MAKEFILE "F77LAGS     = -fast\n";
print MAKEFILE "F77LIBS     = -llapack -lpgplot -lU77 \n";
print MAKEFILE "F77LDFLAGS  =  -L\$(HOME)/Fortran/Lib/\n\n";
#----  F90 -----
print MAKEFILE "F90         = gfortran\n";
print MAKEFILE "F90FLAGS    = -Wall -Wextra\n";
print MAKEFILE "F90LIBS     =  \n";
print MAKEFILE "#F90LIBS     = -lmod -lnr -llapack90 -llapack -lpgplot -lU77 \n";
print MAKEFILE "F90LDFLAGS  = -s -L\$(HOME)/Fortran/Lib/ \n\n";

# Set the LIB and LDFLAGS variable
if (&LanguageCompiler($ARGV[1], @srcs) eq "C") { 
		print MAKEFILE "LIBS        =  \$(CLIBS) \t\n\n"; 
		print MAKEFILE "LDFLAGS     =  \$(CLDFLAGS) \t\n\n"; 
};
if (&LanguageCompiler($ARGV[1], @srcs) eq "CPP") { 
		print MAKEFILE "LIBS        =  \$(CPPLIBS) \t\n\n"; 
		print MAKEFILE "LDFLAGS     =  \$(CPPLDFLAGS) \t\n\n"; 
};
if (&LanguageCompiler($ARGV[1], @srcs) eq "F77"){ 
		print MAKEFILE "LIBS        =  \$(F77LIBS) \t\n\n"; 
		print MAKEFILE "LDFLAGS     =  \$(F77LDFLAGS) \t\n\n"; 
};
if (&LanguageCompiler($ARGV[1], @srcs) eq "F90") { 
		print MAKEFILE "LIBS        =  \$(F90LIBS) \t\n\n"; 
		print MAKEFILE "LDFLAGS     =  \$(F90LDFLAGS) \t\n\n"; 
};


# --- etags ---
print MAKEFILE "ETAGS       = etags\n";

# --- f90 dependencies ---
print MAKEFILE "DEPEND      = makedepf90\n";

# --- svn version ---
print MAKEFILE "SVNVERSION  = svnversion\n\n\n\n";
#print MAKEFILE "VERSION     = VERSION\n\n\n";


#
# make
#
print MAKEFILE "all: etags depend \$(PROG)\n\n";
print MAKEFILE "\$(PROG): \$(OBJS)\n";
print MAKEFILE "\t\$(", &LanguageCompiler($ARGV[1], @srcs);
print MAKEFILE ") \$(LDFLAGS) -o \$@ \$(OBJS) \$(LIBS)\n\n";
#
# make clean
#
print MAKEFILE "clean:\n";
print MAKEFILE "\trm -f \$(PROG) \$(OBJS) *.mod  TAGS VERSION .depend\n\n";
#
# make etags
#
print MAKEFILE "etags:\n";
print MAKEFILE "\t\$(ETAGS) \$(SRCS) *.h90\n\n";
#
# make depend 
#        
print MAKEFILE "depend .depend:\n";
print MAKEFILE "\trm -f .depend\n";
print MAKEFILE "\t\$(DEPEND) \$(SRCS) > .depend \n\n";
#
# make svnversion
#        
print MAKEFILE "version:\n";
print MAKEFILE "\t\$(SVNVERSION) > VERSION \n\n\n";




#
# Make .f90 a valid suffix
#
print MAKEFILE ".SUFFIXES: \$(SUFFIXES) .f90 .F90 .mod\n\n";
#
# .f90 -> .o
#
print MAKEFILE ".f90.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# .f -> .o
#
print MAKEFILE ".f.o:\n";
print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c \$<\n\n";
#
# .c -> .o
#
print MAKEFILE ".c.o:\n";
print MAKEFILE "\t\$(CC) \$(CFLAGS) -c  \$<\n\n";
#
# .cc -> .o
#
print MAKEFILE ".cc.o:\n";
print MAKEFILE "\t\$(CPP) \$(CPPFLAGS) -c \$<\n\n";
#
# .cpp -> .o
#
print MAKEFILE ".cpp.o:\n";
print MAKEFILE "\t\$(CPP) \$(CPPFLAGS) -c \$<\n\n";

print MAKEFILE "# Override the modula2 default built-in rule! \n";
print MAKEFILE "#   (without it, make will try to execute m2c .....) \n";
print MAKEFILE "%.o:%.mod \n";



print MAKEFILE "\n\n\# ... Dependencies\n";
print MAKEFILE "\# ......................\n\n";
print MAKEFILE "\# Include the dependency list (created by \$(DEPEND))\n"; 
print MAKEFILE "include .depend\n\n";







#
# Dependency listings
#
# -----  These routine are currently not in use (IS Oct. 2008) -----
#&MakeDependsf90($ARGV[1]);
#&MakeDepends("*.f *.F",        '^\s*include\s+["\']([^"\']+)["\']');
#&MakeDepends("*.c",            '^\s*#\s*include\s+["\']([^"\']+)["\']');
#&MakeDepends("*.cc *.cpp *.C", '^\s*#\s*include\s+["\']([^"\']+)["\']');


# ---------------------------------------------




#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
#
sub PrintWords {
	 local($columns) = 78 - shift(@_);
	 local($extratab) = shift(@_);
	 local($wordlength);
	 #
	 print MAKEFILE @_[0];
	 $columns -= length(shift(@_));
	 foreach $word (@_) {
			$wordlength = length($word);
			if ($wordlength + 1 < $columns) {
				 print MAKEFILE " $word";
				 $columns -= $wordlength + 1;
				 }
			else {
				 #
				 # Continue onto a new line
				 #
				 if ($extratab) {
						print MAKEFILE " \\\n\t\t$word";
						$columns = 62 - $wordlength;
						}
				 else {
						print MAKEFILE " \\\n\t$word";
						$columns = 70 - $wordlength;
						}
				 }
			}
	 }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
	 local($compiler) = &toLower(shift(@_));
	 local(@srcs) = @_;
	 #
	 if (length($compiler) > 0) {
			CASE: {
				 grep(/^$compiler$/, ("fc", "f77")) &&
						do { $compiler = "F77"; last CASE; };
				 grep(/^$compiler$/, ("c"))   &&
						do { $compiler = "C"; last CASE; };
				 grep(/^$compiler$/, ("cc", "cpp"))   &&
						do { $compiler = "CPP"; last CASE; };
				 $compiler = "F90";
				 }
			}
	 else {
			CASE: {
				 grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
				 grep(/\.(f|F)$/, @srcs) && do { $compiler = "F77";  last CASE; };
				 grep(/\.c$/, @srcs)     && do { $compiler = "C";  last CASE; };
		 grep(/\.cc$/, @srcs)    && do { $compiler = "CPP"; last CASE; };
		 grep(/\.cpp$/, @srcs)   && do { $compiler = "CPP"; last CASE; };
		 grep(/\.C$/, @srcs)     && do { $compiler = "CPP"; last CASE; };
				 $compiler = "???";
				 }
			}
	 $compiler;
	 }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
	 local($string) = @_[0];
	 $string =~ tr/A-Z/a-z/;
	 $string;
	 }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
	 local(@words);
	 foreach $word (@_) {
			if ($word ne $words[$#words]) {
				 push(@words, $word);
				 }
			}
	 @words;
	 }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
	 local(@incs);
	 local($lang) = @_[0];
	 local($pattern) = @_[1];
	 #
	 foreach $file (<${lang}>) {
			open(FILE, $file) || warn "Cannot open $file: $!\n";
			while (<FILE>) {
				 /$pattern/i && push(@incs, $1);
				 }
			if (defined @incs) {
				 $file =~ s/\.[^.]+$/.o/;
				 print MAKEFILE "$file: ";
				 &PrintWords(length($file) + 2, 0, @incs);
				 print MAKEFILE "\n";
				 undef @incs;
				 }
			}
	 }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
	 local($compiler) = &toLower(@_[0]);
	 local(@dependencies);
	 local(%filename);
	 local(@incs);
	 local(@modules);
	 local($objfile);
	 #
	 # Associate each module with the name of the file that contains it
	 #
	 foreach $file (<*.f90>) {
			open(FILE, $file) || warn "Cannot open $file: $!\n";
			while (<FILE>) {
				 /^\s*module\s+([^\s!]+)/i &&
						($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
				 }
			}
	 #
	 # Print the dependencies of each file that has one or more include's or
	 # references one or more modules
	 #
	 foreach $file (<*.f90>) {
			open(FILE, $file);
			while (<FILE>) {
				 /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
				 /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
				 }
			if (defined @incs || defined @modules) {
				 ($objfile = $file) =~ s/\.f90$/.o/;
				 print MAKEFILE "$objfile: ";
				 undef @dependencies;
				 foreach $module (@modules) {
						push(@dependencies, $filename{$module});
						}
				 @dependencies = &uniq(sort(@dependencies));
				 &PrintWords(length($objfile) + 2, 0,
										 @dependencies, &uniq(sort(@incs)));
				 print MAKEFILE "\n";
				 undef @incs;
				 undef @modules;
				 #
				 # Cray F90 compiler
				 #
				 if ($compiler eq "cray") {
						print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c -O";
						foreach $depend (@dependencies) {
							 push(@modules, "-p", $depend);
							 }
						push(@modules, $file);
						&PrintWords(30, 1, @modules);
						print MAKEFILE "\n";
						undef @modules;
						}
				 #
				 # ParaSoft F90 compiler
				 #
				 if ($compiler eq "parasoft") {
						print MAKEFILE "\t\$(F90) \$(F90FLAGS) -c -O";
						foreach $depend (@dependencies) {
							 $depend =~ s/\.o$/.f90/;
							 push(@modules, "-module", $depend);
							 }
						push(@modules, $file);
						&PrintWords(30, 1, @modules);
						print MAKEFILE "\n";
						undef @modules;
						}
				 }
			}
	 }