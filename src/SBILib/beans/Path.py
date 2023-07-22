"""Path

author: jbonet
date:   02/2013

@oliva's lab
"""

"""
Import Standard Libraries
"""
import os
import fnmatch
import shutil

"""
Dependences in SBI library
"""
from SBILib.error import FileError
from .             import File
from SBILib           import SBIglobals

class Path(object):
    """Path

        Path is a collection of @staticmethod (can be called without declaring an instance of the class) to 
        help in file and directory management.

        Methods:
            > list_files()      : Returns any file in a directory tree matching a specific pattern.
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - pattern           (string): Expression to match (ls-like format)
                                                                  (Accepts list of strings)
                                                                  @Default: *
                                    - avoid_empty_files (bool)  : Ignore files with size 0
                                                                  @Default: True
                                    - rootless          (bool)  : When False, the name of the files are returned with
                                                                  absolute path. Otherwise the root is removed.
                                                                  @Default: False
                                  @Yields file names

            > list_directories(): Returns all dirs in a directory tree.
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - rootless          (bool)  : When False, the name of the dirs are returned with
                                                                  absolute path. Otherwise the root is removed.
                                                                  @Default: False
                                  @Yields directory names

            > sync_directories(): Syncronizes two directory trees.
                                    - sourcedir      (string): Name of the original directory
                                                               @Mandatory
                                    - destinationdir (string): Name of the syncronized directory
                                                               @Mandatory
                                    - by_dir         (bool)  : Sync is performed at directory level
                                                               @Default: True
                                    - by_file        (bool)  : Sync is performed at file level
                                                               @Default: False
                                    - listdif        (bool)  : When True, sync is not performed, only listed
                                                               @Default: False
                                    @Yields file/dir names when listdif is True
                                    @Raises AttributeError if by_dir and by_file are both True

            > do_files_match()  : Determines if there is any file in a directory tree matching a specific pattern.
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - pattern           (string): Expression to match (ls-like format)
                                                                  @Default: *
                                    - avoid_empty_files (bool)  : Ignore files with size 0
                                                                  @Default: True
                                  @Returns Boolean

            - mkdir()           : Creates a new directory. 
                                  Ignores the command if the directory already exists.
                                  Creates parent directories as required.
                                  Recipe: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/82465
                                    - newdir (string): Name of the new directory.
                                                       @Mandatory
                                  @Raises FileError

            - copy_file()       : Copy a file from one place to another.
                                    - source      (string): Name of the original file
                                                            @Mandatory
                                    - destination (string): Name of the new copied file
                                                            @Mandatory
                                    -overwrite    (bool)  : Ignores if a file with the same name existed before
                                                            @Default: False
                                  @Raises FileError
    """
    @staticmethod
    def list_files(root = os.curdir, pattern = '*', avoid_empty_files = True, rootless = False):
        """
            > list_files()      : Returns any file in a directory tree matching a specific pattern
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - pattern           (string): Expression to match (ls-like format)
                                                                  (Accepts list of strings)
                                                                  @Default: *
                                    - avoid_empty_files (bool)  : Ignore files with size 0
                                                                  @Default: True
                                    - rootless          (bool)  : When False, the name of the files are returned with
                                                                  absolute path. Otherwise the root is removed.
                                                                  @Default: False
                                  @Yields file names
        """
        if os.path.isfile(root): yield root

        search_patterns = []
        if not isinstance(pattern, list):
            search_patterns.append(pattern)
        else:
            search_patterns = pattern
        for pat in search_patterns:
            for path, dirs, files in os.walk(os.path.abspath(root)):
                for filename in fnmatch.filter(files, pat):
                    if not avoid_empty_files or os.path.getsize(os.path.join(path, filename)) > 0:
                        SBIglobals.alert('debug', Path(), 'Found file {0}'.format(os.path.join(path,filename)))
                        if not rootless:
                            yield os.path.join(path, filename)
                        else:
                            root = os.path.abspath(root) + "/"
                            yield os.path.join(path, filename).replace(root,'')

    @staticmethod
    def list_directories(root = os.curdir, rootless = False):
        """
            > list_directories(): Returns all dirs in a directory tree
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - rootless          (bool)  : When False, the name of the dirs are returned with
                                                                  absolute path. Otherwise the root is removed.
                                                                  @Default: False
                                  @Yields directory names
        """
        for path, dirs, files in os.walk(os.path.abspath(root)):
            for onedir in dirs:
                SBIglobals.alert('debug', Path(), 'Found directory {0}'.format(os.path.join(path,onedir)))
                if not rootless:
                    yield os.path.join(path,onedir)
                else:
                    yield os.path.join(path,onedir).replace(root, '')

    @staticmethod
    def sync_directories(sourcedir, destinationdir, by_dir = True, by_file = False, listdif = False):
        """
            > sync_directories(): Syncronizes two directory trees.
                                    - sourcedir      (string): Name of the original directory
                                                               @Mandatory
                                    - destinationdir (string): Name of the syncronized directory
                                                               @Mandatory
                                    - by_dir         (bool)  : Sync is performed at directory level
                                                               @Default: True
                                    - by_file        (bool)  : Sync is performed at file level
                                                               @Default: False
                                    - listdif        (bool)  : When True, sync is not performed, only listed
                                                               @Default: False
                                    @Yields file/dir names when listdif is True
                                    @Raises AttributeError if by_dir and by_file are both True
        """
        if by_file is True and by_dir is True:
            raise AttributeError('Both sync methods can not be active simultaneously')

        source_dirs = set()
        if by_dir:
            for onedir in Path.list_directories(root = sourcedir, rootless = True):
                source_dirs.add(onedir)
            for onedir in Path.list_directories(root = destinationdir, rootless = True):
                if onedir in source_dirs:
                    source_dirs.remove(onedir)
        if by_file:
            for onedir in Path.list_files(root = sourcedir, rootless = True):
                source_dirs.add(onedir)
            for onedir in Path.list_files(root = destinationdir, rootless = True):
                if onedir in source_dirs:
                    source_dirs.remove(onedir)

        previous_dir   = '#not_a_dir_at_all'
        sourcedir      = os.path.abspath(sourcedir)
        destinationdir = os.path.abspath(destinationdir)
        for onedir in sorted(source_dirs):
            onedir = onedir.lstrip('/')
            if not os.path.commonprefix([previous_dir, onedir]) == previous_dir:
                SBIglobals.alert('verbose', Path(), '{0} is different from {1} to {2}'.format(onedir,sourcedir,destinationdir))
                if not listdif:
                    fullori = os.path.join(sourcedir,onedir)
                    shutil.copytree(os.path.join(sourcedir,onedir),os.path.join(destinationdir,onedir))
                else:
                    yield os.path.join(sourcedir, onedir)
                previous_dir = onedir

    @staticmethod
    def do_files_match(root = os.curdir, pattern = '*', avoid_empty_files = True):
        """
            > do_files_match()  : Determines if there is any file in a directory tree matching a specific pattern
                                    - root              (string): Root of the directory tree to search
                                                                  @Default: current working directory
                                    - pattern           (string): Expression to match (ls-like format)
                                                                  @Default: *
                                    - avoid_empty_files (bool)  : Ignore files with size 0
                                                                  @Default: True
                                  @Returns Boolean
        """
        if os.path.isfile(root): return True

        search_patterns = []
        if not isinstance(pattern, list):
            search_patterns.append(pattern)
        else:
            search_patterns = pattern
        for pat in search_patterns:
            for path, dirs, files in os.walk(os.path.abspath(root)):
                for filename in fnmatch.filter(files, pat):
                    if not avoid_empty_files or os.path.getsize(os.path.join(path, filename)) > 0:
                        return True
        return False

    @staticmethod
    def mkdir(newdir):
        """
            - mkdir()           : Creates a new directory. 
                                  Ignores the command if the directory already exists.
                                  Creates parent directories as required.
                                  Recipe: http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/82465
                                    - newdir (string): Name of the new directory.
                                                       @Mandatory
                                  @Raises FileError
        """
        if os.path.isdir(newdir): 
            SBIglobals.alert('debug', Path(), 'A directory with that name already exists.')
            pass
        elif os.path.isfile(newdir):
            raise FileError(3,newdir,'filewithname');
        else:
            SBIglobals.alert('debug', Path(), 'Creating dir {0}.'.format(newdir))
            head, tail = os.path.split(newdir)
            if head and not os.path.isdir(head):
                Path.mkdir(head)
            if tail:
                os.mkdir(newdir)

    @staticmethod
    def copy_file(source, destination, overwrite = False):
        """
            - copy_file()       : Copy a file from one place to another.
                                    - source      (string): Name of the original file
                                                            @Mandatory
                                    - destination (string): Name of the new copied file
                                                            @Mandatory
                                    -overwrite    (bool)  : Ignores if a file with the same name existed before
                                                            @Default: False
                                  @Raises FileError
        """
        # Local overwrite takes precedence over Global overwrite
        overwrite = SBIglobals.decide_overwrite(overwrite)

        if not os.path.exists(source):
            raise FileError(3, source, 'noexists')
        if os.path.exists(destination) and not overwrite:
            raise FileError(3, destination, 'exists')

        shutil.copyfile(src = source, dst = destination)
