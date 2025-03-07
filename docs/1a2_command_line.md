---
layout: default
title: The command line
nav_order: 5
parent: 1. Introduction
description: Introduction to the command line
published: true
---

FINAL VERSION
{: .label .label-green }


This is an introduction to the *UNIX* command line.
It covers some commands and concepts that are widely used in the UNIX command line.
Examples are provided throughout.

<br>
<details open markdown="block">
  <summary>
    <strong>Table of contents</strong>
  </summary>
  {: .text-delta }
- TOC
{:toc}
</details>

<!--
## Contents

- [Opening the terminal](#opening-the-terminal)
- [The working directory](#the-working-directory)
- [Cloning the `bioinfo-notebook` project into your home directory](#cloning-the-bioinfo-notebook-project-into-your-home-directory)
- [Changing working directories](#changing-working-directories)
- [Comments and broken lines](#comments-and-broken-lines)
- [Listing directory content with the `ls` command](#listing-directory-content-with-the-ls-command)
- [Relative paths](#relative-paths)
- [Making directories with `mkdir`](#making-directories-with-mkdir)
- [Removing files and directories with the `rm` command](#removing-files-and-directories-with-the-rm-command)
- [Using the `head` command](#using-the-head-command)
- [The `tail` command](#the-tail-command)
- [Using the `--help` argument](#using-the---help-argument)
- [Downloading files with `wget`](#downloading-files-with-wget)
- [Moving and copying with `mv` and `cp`](#moving-and-copying-with-mv-and-cp)
- [Running the Linux setup shell script](#running-the-linux-setup-shell-script)
    - [Video demonstration](#video-demonstration)
- [Exercise](#exercise)
- [See also](#see-also)
-->

## The terminal

The terminal is the window in which the command line runs.

The UNIX command line will look like this:

```
(UNIX username)@(Computer's alias):~$ _
```

This is called the *bash prompt*...

- *UNIX username* is the your system username.
- *Computer's alias* is the name UNIX uses to refer to your computer. This will likely contain the model of your computer (e.g. `Latitude-E7270`).
- The tilde (`~`) indicates that your home directory is the current *working directory*. The home directory is located at `/home/` followed by your UNIX username.
- The dollar sign (`$`) indicates that the terminal is using the `bash` shell language.

In this page, `st24_16@banjo:~$` will be used as an example bash prompt.

## The working directory

A *directory* is the same as a folder on your desktop.
If you have a "Pictures" folder on your computer's desktop, this folder is a directory within the desktop directory: its [path](#relative-paths) is `Desktop/Pictures/`.
The *working directory* is the directory that the command line is currently using.
Any files that you create will be in this directory, and any commands you run will use files in this directory (unless you use a [path](#relative-paths)).

To see the current working directory, type `pwd` into the command line and press `Enter` (or `Return`).
This will run the "print working directory" command, which will print the path to the current directory in the terminal.
In the context of command line programs, "print" means "display in the terminal".

```bash
st24_16@banjo:~$ pwd
```

This will give:

```bash
st24_16@banjo:~$ pwd
/home/st24_16
```


## Changing working directories

The working directory of the command line can be changed using the `cd` (change directory) command.
Use `pwd` before and after the `cd` command to check which directory you have moved from and to.
Type `cd /data2` to change the working directory to the `/data2` directory.

```bash
st24_16@banjo:~$ pwd
/home/st24_16
st24_16@banjo:~$ cd /data2
st24_16@banjo:/data2$ pwd
/data2
```

You do not need to type this command out in full, typing `cd /d` and pressing the `Tab` key should complete this command.
Notice that the working directory part of the bash prompt has changed from `~` to `~/data2`.

When you press `Tab` in the command line, it will try to complete the file or directory name for you.
If there are multiple options, these options will be printed in the terminal.
If there is only one possible option, this option will be filled in.
For example, pressing `Tab` twice after typing `cd /d` could give the following directories beginning with "d" as options:

```bash
st24_16@banjo:~$ cd /d
data/  data2/ dev/
```

In UNIX, directory names are case-sensitive.
This means that `~/downloads/` is a different directory than `~/Downloads/`.

The command `cd ../` can be used to move to move up one directory, `cd ../../` can be used to move up two directories, etc.

In this example, moving up one directory from `/data2/biotecnologie_molecolari_magris/epigenomics/` changes the working directory to `/data2/biotecnologie_molecolari_magris/`.

Moving up two directories from `/data2/biotecnologie_molecolari_magris/epigenomics/` changes the working directory to `/data2`.

```bash
cd /data2/biotecnologie_molecolari_magris/epigenomics/
st24_16@banjo:/data2/biotecnologie_molecolari_magris/epigenomics$ pwd
/data2/biotecnologie_molecolari_magris/epigenomics

# go up one directory
cd ../

# get the present working directory 
pwd
/data2/biotecnologie_molecolari_magris

# go down one directory
cd epigenomics/

# go up two directories
cd ../../

# get the present working directory 
pwd
/data2/
```

## Listing directory content with the `ls` command

The `ls` command can be used to list the files and directories within the current working directory.
Try using the `ls` command in the `/data2/` directory.

```bash
ls
biotecnologie_molecolari_depaoli  biotecnologie_molecolari_magris
```

## Relative paths

From the `/data2/student_space/st24_16_folder/` working directory, all of the files used during the tutorials can be accessed using *relative paths*.
Without paths, the program only knows which file or directory to look for **in the working directory.**
With paths, the program knows which file or directory to look for *and* how to get to it from the working directory.

The `ls` command can be used to list the content of the folder `/data2/biotecnologie_molecolari_magris/` from the `/data2` working directory using `ls biotecnologie_molecolari_magris/`.

```bash
pwd
/data2
```
```bash
ls biotecnologie_molecolari_magris/
epigenomics
```

## Using the `--help` argument

For a lot of command line programs, using the command with the `--help` (or `-h`) argument will display a short help message on how that command is used.

```bash
ls --help
Usage: ls [OPTION]... [FILE]...
List information about the FILEs (the current directory by default).
Sort entries alphabetically if none of -cftuvSUX nor --sort is specified.

Mandatory arguments to long options are mandatory for short options too.
  -a, --all                  do not ignore entries starting with .
  -A, --almost-all           do not list implied . and ..
      --author               with -l, print the author of each file
  -b, --escape               print C-style escapes for nongraphic characters
      --block-size=SIZE      scale sizes by SIZE before printing them; e.g.,
                               '--block-size=M' prints sizes in units of
                               1,048,576 bytes; see SIZE format below
  -B, --ignore-backups       do not list implied entries ending with ~
  -c                         with -lt: sort by, and show, ctime (time of last
                               modification of file status information);
                               with -l: show ctime and sort by name;
                               otherwise: sort by ctime, newest first
  -C                         list entries by columns
      --color[=WHEN]         colorize the output; WHEN can be 'never', 'auto',
                               or 'always' (the default); more info below
  -d, --directory            list directories themselves, not their contents
  -D, --dired                generate output designed for Emacs' dired mode
  -f                         do not sort, enable -aU, disable -ls --color
  -F, --classify             append indicator (one of */=>@|) to entries
      --file-type            likewise, except do not append '*'
      --format=WORD          across -x, commas -m, horizontal -x, long -l,
                               single-column -1, verbose -l, vertical -C
      --full-time            like -l --time-style=full-iso
  -g                         like -l, but do not list owner
      --group-directories-first
                             group directories before files;
                               can be augmented with a --sort option, but any
                               use of --sort=none (-U) disables grouping
  -G, --no-group             in a long listing, don't print group names
  -h, --human-readable       with -l and/or -s, print human readable sizes
                               (e.g., 1K 234M 2G)
      --si                   likewise, but use powers of 1000 not 1024
  -H, --dereference-command-line
                             follow symbolic links listed on the command line
      --dereference-command-line-symlink-to-dir
                             follow each command line symbolic link
                               that points to a directory
      --hide=PATTERN         do not list implied entries matching shell PATTERN
                               (overridden by -a or -A)
      --indicator-style=WORD  append indicator with style WORD to entry names:
                               none (default), slash (-p),
                               file-type (--file-type), classify (-F)
  -i, --inode                print the index number of each file
  -I, --ignore=PATTERN       do not list implied entries matching shell PATTERN
  -k, --kibibytes            default to 1024-byte blocks for disk usage
  -l                         use a long listing format
  -L, --dereference          when showing file information for a symbolic
                               link, show information for the file the link
                               references rather than for the link itself
  -m                         fill width with a comma separated list of entries
  -n, --numeric-uid-gid      like -l, but list numeric user and group IDs
  -N, --literal              print raw entry names (don't treat e.g. control
                               characters specially)
  -o                         like -l, but do not list group information
  -p, --indicator-style=slash
                             append / indicator to directories
  -q, --hide-control-chars   print ? instead of nongraphic characters
      --show-control-chars   show nongraphic characters as-is (the default,
                               unless program is 'ls' and output is a terminal)
  -Q, --quote-name           enclose entry names in double quotes
      --quoting-style=WORD   use quoting style WORD for entry names:
                               literal, locale, shell, shell-always, c, escape
  -r, --reverse              reverse order while sorting
  -R, --recursive            list subdirectories recursively
  -s, --size                 print the allocated size of each file, in blocks
  -S                         sort by file size
      --sort=WORD            sort by WORD instead of name: none (-U), size (-S),
                               time (-t), version (-v), extension (-X)
      --time=WORD            with -l, show time as WORD instead of default
                               modification time: atime or access or use (-u)
                               ctime or status (-c); also use specified time
                               as sort key if --sort=time
      --time-style=STYLE     with -l, show times using style STYLE:
                               full-iso, long-iso, iso, locale, or +FORMAT;
                               FORMAT is interpreted like in 'date'; if FORMAT
                               is FORMAT1<newline>FORMAT2, then FORMAT1 applies
                               to non-recent files and FORMAT2 to recent files;
                               if STYLE is prefixed with 'posix-', STYLE
                               takes effect only outside the POSIX locale
  -t                         sort by modification time, newest first
  -T, --tabsize=COLS         assume tab stops at each COLS instead of 8
  -u                         with -lt: sort by, and show, access time;
                               with -l: show access time and sort by name;
                               otherwise: sort by access time
  -U                         do not sort; list entries in directory order
  -v                         natural sort of (version) numbers within text
  -w, --width=COLS           assume screen width instead of current value
  -x                         list entries by lines instead of by columns
  -X                         sort alphabetically by entry extension
  -Z, --context              print any security context of each file
  -1                         list one file per line
      --help     display this help and exit
      --version  output version information and exit

The SIZE argument is an integer and optional unit (example: 10K is 10*1024).
Units are K,M,G,T,P,E,Z,Y (powers of 1024) or KB,MB,... (powers of 1000).

Using color to distinguish file types is disabled both by default and
with --color=never.  With --color=auto, ls emits color codes only when
standard output is connected to a terminal.  The LS_COLORS environment
variable can change the settings.  Use the dircolors command to set it.

Exit status:
 0  if OK,
 1  if minor problems (e.g., cannot access subdirectory),
 2  if serious trouble (e.g., cannot access command-line argument).

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/ls>
or available locally via: info '(coreutils) ls invocation'
g
```

For more in-depth help with a command, the `man` (manual) command is useful, if it is available.
This will open a manual for the program within the terminal.
The `Up` and `Down` arrow keys (or `j` and `k`) can be used to scroll through these manuals, and `q` is used to exit.

```bash
man ls

LS(1)                        User Commands                                  LS(1)

NAME
       ls - list directory contents

SYNOPSIS
       ls [OPTION]... [FILE]...

DESCRIPTION
       List information about the FILEs (the current directory by default).  Sort entries alphabetically if none of -cftuvSUX nor --sort is spec‐
       ified.

```


## Creating directories with `mkdir`

The command `mkdir` can be used to create a directory.

```bash
mkdir --help
Usage: mkdir [OPTION]... DIRECTORY...
Create the DIRECTORY(ies), if they do not already exist.

Mandatory arguments to long options are mandatory for short options too.
  -m, --mode=MODE   set file mode (as in chmod), not a=rwx - umask
  -p, --parents     no error if existing, make parent directories as needed
  -v, --verbose     print a message for each created directory
  -Z                   set SELinux security context of each created directory
                         to the default type
      --context[=CTX]  like -Z, or if CTX is specified then set the SELinux
                         or SMACK security context to CTX
      --help     display this help and exit
      --version  output version information and exit

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/mkdir>
or available locally via: info '(coreutils) mkdir invocation'
```

In this example, `mkdir` is used to create a directory called `example`, which is located in (and is therefore a subdirectory of) the working directory (`/data2/student_space/st24_16_folder`).
This is the equivalent of making a new folder within a folder on your desktop.

```bash
cd /data2/student_space/st24_16_folder/

mkdir example

cd example/

pwd

/data2/student_space/st24_16_folder/example
```

## Removing files and directories with the `rm` command

The `rm` (remove) command is used to remove files and directories.
This command must be used *with care,* as it is not the same as deleting a file or directory on a desktop.
Deleted files on a desktop are moved to a Recycle Bin or Trash folder, and are only permanently deleted once they are removed from this folder.
When using the `rm` command, files and directories are permanently deleted; they cannot be recovered.

Here is the help text for the `rm` command:

```bash
rm --help
Usage: rm [OPTION]... [FILE]...
Remove (unlink) the FILE(s).

  -f, --force           ignore nonexistent files and arguments, never prompt
  -i                    prompt before every removal
  -I                    prompt once before removing more than three files, or
                          when removing recursively; less intrusive than -i,
                          while still giving protection against most mistakes
      --interactive[=WHEN]  prompt according to WHEN: never, once (-I), or
                          always (-i); without WHEN, prompt always
      --one-file-system  when removing a hierarchy recursively, skip any
                          directory that is on a file system different from
                          that of the corresponding command line argument
      --no-preserve-root  do not treat '/' specially
      --preserve-root   do not remove '/' (default)
  -r, -R, --recursive   remove directories and their contents recursively
  -d, --dir             remove empty directories
  -v, --verbose         explain what is being done
      --help     display this help and exit
      --version  output version information and exit

By default, rm does not remove directories.  Use the --recursive (-r or -R)
option to remove each listed directory, too, along with all of its contents.

To remove a file whose name starts with a '-', for example '-foo',
use one of these commands:
  rm -- -foo

  rm ./-foo

Note that if you use rm to remove a file, it might be possible to recover
some of its contents, given sufficient expertise and/or time.  For greater
assurance that the contents are truly unrecoverable, consider using shred.

GNU coreutils online help: <http://www.gnu.org/software/coreutils/>
Full documentation at: <http://www.gnu.org/software/coreutils/rm>
or available locally via: info '(coreutils) rm invocation'
```

To remove files, simply type `rm` followed by the name of the file, and press `Enter`/`Return`:

```bash
pwd
/data2/student_space/st24_16_folder/example

# Create an empty file using touch
touch example_file.txt

# List the content of the current directory 
ls
example_file.txt

# Remove the file 
rm example_file.txt 

# List again the content of the current directory 
ls
```

To remove directories, `rm -r` must be typed, followed by the name of the directory.
**Use the `rm -r` command with caution**, as it will permanently delete entire directories and *all* of their contents.
This includes files and other directories within this directory.

In this example, `rm -r` is used to remove the `example/` directory from the home directory:

```bash
cd ../

pwd

rm -r example/

cd example

bash: cd: example: No such file or directory

```

## Using the `head` command

The `head` command can be used to view the first part (the head) of a file or files.
This command is useful for examining very large files quickly.
From the current directory, use `head /data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta` to view the head of that FASTA file.

```bash
head /data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta
>CASBKX010000001.1 Saccharomyces cerevisiae genome assembly, contig: chrI, whole genome shotgun sequence
CACACACACCACACCCACACCCACACACCCACACCCACACACCCACACCACACCCACACCACACCCACAC
CACACCACACCCACACCACACACCACACCCACACACCACACCACACCCACACCACACCCACACCACACCC
ACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACCACACCCAC
ACACCACACCCACACCCACACACACACACACCACACCCACACCCACACCCACACACCCACACACCCACAA
CCACACCCACACCACACCCACACCCACACACACCACACACCACACCCACACACACCCACACAACCCCGCT
```

The `head` command takes optional *arguments*.
Command line arguments are extra instructions given to a program when it runs.
One of the arguments that can be given to `head` is `-n`, which specifies how many lines we want to print per file.
The command `head -n 5 data/example_nucleotide_sequence.fasta` prints the first 5 lines of the file `example_nucleotide_sequence.fasta`.

```bash
head -n 5 /data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta
>CASBKX010000001.1 Saccharomyces cerevisiae genome assembly, contig: chrI, whole genome shotgun sequence
CACACACACCACACCCACACCCACACACCCACACCCACACACCCACACCACACCCACACCACACCCACAC
CACACCACACCCACACCACACACCACACCCACACACCACACCACACCCACACCACACCCACACCACACCC
ACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACCACACCCAC
ACACCACACCCACACCCACACACACACACACCACACCCACACCCACACCCACACACCCACACACCCACAA

```

## The `tail` command

The `tail` command is the equivalent of the command `head`, but prints the last part of files instead.
This is useful for quickly examining the information at the ends of files.

```bash
tail -n 5 /data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta
CACACACACCACACCCACACCCACACACCCACACCCACACACCCACACCACACCCACACCACACCCACAC
CACACCACACCCACACCACACACCACACCCACACACCACACCACACCCACACCACACCCACACCACACCC
ACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACCACACCCAC
ACACCACACCCACACCCACACACACACACACCACACCCACACCCACACCCACACACCCACACACCCACAA
CCACACCCACACCACACCCACACCCACACACACCACACACCACACCCACACACACCCACACAACCCCGCT
```

## Comments and broken lines

The number sign or hash (`#`) is used for *comments* in bash.
Comments are used in programming to explain or annotate code; they are not run as code.
Anything written after a `#` symbol in the command line is a comment.

```bash
pwd # This prints the working directory
/data2/
```

Long lines of code can be hard to read.
For legibility, long lines of command line code can be broken up with the backslash (`\`).
The `\` characters split lines to make them more readable, they do not affect how the code functions.

```bash
head /data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta

>CASBKX010000001.1 Saccharomyces cerevisiae genome assembly, contig: chrI, whole genome shotgun sequence
CACACACACCACACCCACACCCACACACCCACACCCACACACCCACACCACACCCACACCACACCCACAC
CACACCACACCCACACCACACACCACACCCACACACCACACCACACCCACACCACACCCACACCACACCC
ACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACCACACCCAC
ACACCACACCCACACCCACACACACACACACCACACCCACACCCACACCCACACACCCACACACCCACAA
CCACACCCACACCACACCCACACCCACACACACCACACACCACACCCACACACACCCACACAACCCCGCT

# The same can be achieved 
head \
/data2/biotecnologie_molecolari_magris/epigenomics/bash/example_sequence.fasta

>CASBKX010000001.1 Saccharomyces cerevisiae genome assembly, contig: chrI, whole genome shotgun sequence
CACACACACCACACCCACACCCACACACCCACACCCACACACCCACACCACACCCACACCACACCCACAC
CACACCACACCCACACCACACACCACACCCACACACCACACCACACCCACACCACACCCACACCACACCC
ACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACCACACCCAC
ACACCACACCCACACCCACACACACACACACCACACCCACACCCACACCCACACACCCACACACCCACAA
CCACACCCACACCACACCCACACCCACACACACCACACACCACACCCACACACACCCACACAACCCCGCT
```


## Downloading files with `wget`

The `wget` command can be used to download files from the directory using the command line.
It can take a URL as an argument, and download the file found at that URL into the working directory.

In this example, `wget` is used to download the amino acid sequence of [haemoglobin subunit alpha from UniProt](https://www.uniprot.org/uniprot/P69905):

```bash
cd /data2/student_space/st24_16_folder/

# Create the directory 
mkdir -p epigenomics/bash

# Move to the directory of interest 
cd epigenomics/bash

# Download the insulin protein sequence 
wget \
https://www.uniprot.org/uniprot/P01308.fasta
--2025-02-03 16:28:20--  https://www.uniprot.org/uniprot/P01308.fasta
Resolving www.uniprot.org (www.uniprot.org)... 193.62.193.81
Connecting to www.uniprot.org (www.uniprot.org)|193.62.193.81|:443... connected.
HTTP request sent, awaiting response... 301 Moved Permanently

....

P01308.fasta        100%[================================================>]     182  --.-KB/s   in 0s
 

2025-02-03 16:28:21 (5.00 MB/s) - ‘P01308.fasta’ saved [182/182]

```

```bash
head P01308.fasta # show the first rows of the downloaded file
>sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED
LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN
```



<!--
It may be necessary to download sequences also from public repositories (for example fastq files)

# ADD the step of downloading data from ncbi archive -> need to add the powerpoint presentation 
-->

## Moving and copying with `mv` and `cp`

The `mv` and `cp` commands can be used to move and copy files or directories from the command line.
Both of these commands take a target file/directory and a destination file/directory as arguments.
The `cp` commands creates a copy of a file/directory, whereas `mv` moves it to a new destination.
The `mv` command can be used to rename files, by moving it to the same directory, but with a new file name.

```bash
# Create a subfolder
mkdir test 

# Move the file to new folder
mv P01308.fasta test/

# Check if the file has been moved 
ls test/

# Rename the file 
mv test/P01308.fasta test/ins.fasta

# Check if the file has been renamed 
ls test/

# Create a copy of the fasta file in the working directory
cp test/ins.fasta . 

# Check the content of the working directory
ls 

ins.fasta test
# Create a copy of the fasta file in the working directory with a different name
cp test/ins.fasta insulin.fasta

# Check the content of the working directory
ls 

ins.fasta insulin.fasta  test
```

## Compressed files 
It may happen that files downloaded or shared are in a compressed format, such as `.zip`, `.tar`, `.gz`, etc. The most widely used compression format can be decompressed using the following commands:

```bash
gunzip # for .gz files

bunzip2 # for .bz2 files

unzip # for .zip files
```

Viceversa, to compress a file, you can use the following commands:
```bash
gzip # to create .gz files

# or
gzip -c # to create .gz files with a different name (without overwriting the original file)

bzip2 # to create .bz2 files
```

In order to create or manipulate archive of files, you can use the following commands:

```bash
tar -cf archive_name.tar file1 file2 file3 # to create an archive of files named archive_name.tar

tar -zcf archive_name.tar.gz file1 file2 file3 # to create an archive of files named archive_name.tar.gz (compressed format)
```

In order to extract the content of an archive, you can use the following commands:

```bash
tar -xf archive_name.tar # to extract the content of the archive named archive_name.tar


tar -zxf archive_name.tar.gz # to extract the content of the archive named archive_name.tar.gz (compressed format)
```

## Wildcards

In the UNIX terminal, an asterisk (`*`) acts as a *wildcard*.
This means that any files or directories that can replace this character will replace it.

The asterisk command is especially useful for selecting files with the same *file extension*.
The file extension is the part of the filename after the full stop that specifies the file type: for example, a file ending in `.txt` is a text file.
In a directory with many text files, `*.txt` selects all of the text files.


## AWK command line
The basic syntax of the AWK command is: 

```bash
awk 'BEGIN { begin command(s) } { line command(s) } END { end command(s) }' input_file
```

Without BEGIN and END 
```bash
awk '{ line command(s) } ' input_file

# For example

awk '{print $0}' input_file # prints the entire line 
awk '{print $3}' input_file # prints the third column of the line
awk '{print $2"\t"$3}' input_file # prints the second and third columns of the line
# the same can be achieved 
awk 'OFS="\t" {print $2,$3}' input_file # prints the second and third columns of the line with tab as the output field separator
```

Field separator can be indicated in the BEGIN section for the input file (FS), and for the output file field separator (OFS). The column of the file can be referred to by its number (1, 2, 3, etc.) or by its name (if the file has a header) and separated using the comma:

```bash
awk 'BEGIN {FS=","; OFS="\t"} {print $2,$3}' input_file
```

`If` functions can be introduced in the awk command line.
For exaple if we want to select only rows that have a certain value in the third column, we can use the following command:

```bash
awk '{if($3=="chr3") print $0}' input_file

# The same can be achieved
awk '$3=="chr3"{print $0}' input_file

# or 

awk '{if($3=="chr3" && $4==10) print $0}' input_file
```

Multiple `if` commands can be combined
```bash
awk '{if(condition 1) actionA; if(condition 2 && condition 3) actionB}' input_file
```

Let's try with an exercise. Suppose to have a file with the following content:
```bash
Carlo   Biologia        80
Carlo   Chimica 100
Donato  Biologia        100
Donato  Chimica 90
Francesca       Biologia        90
Francesca       Chimica 80
Chiara  Biologia        100
Chiara  Chimica 70
```

File is located here:

`/data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt`

- Discard the names and extract only the 2nd and 3rd columns:

```bash
awk '{print $2"\t"$3}' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt

# The same can be achieved 
awk 'OFS="\t" {print $2,$3}' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt
```

- Extract only the rows with `Biologia`

```bash
awk '{if($2=="Biologia") print $1"\t"$3}' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt

# In a similar manner we can write the code
awk '/Biologia/' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt

# Or

awk '/Biologia/ {print $1"\t"$3}' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt
```

- Extract only the names of student that achieved more or equal to 80 in `Chimica`

```bash
awk '{if($2=="Chimica" && $3>=80) print $1}' /data2/biotecnologie_molecolari_magris/epigenomics/bash/awk_example.txt
```