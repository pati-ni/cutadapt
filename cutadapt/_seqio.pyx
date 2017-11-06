# kate: syntax Python;
# cython: profile=False, emit_code_comments=False
from __future__ import print_function, division, absolute_import
from xopen import xopen
from .seqio import _shorten, FormatError, SequenceReader

from libc.string cimport strncmp

# TODO/ideas
# - Sequence:
#   -non-ASCII qualities could be created (and cached) on access
#   - remove second_header
#   - remove match attribute (?)
# -

cdef class Sequence(object):
	"""
	A record in a FASTA or FASTQ file. For FASTA, the qualities attribute
	is None. For FASTQ, qualities is a string and it contains the qualities
	encoded as ascii(qual+33).

	If an adapter has been matched to the sequence, the 'match' attribute is
	set to the corresponding Match instance.
	"""
	cdef:
		public str name
		public str sequence
		public str qualities
		public bint second_header
		public object match

	def __init__(self, str name, str sequence, str qualities=None, bint second_header=False,
	        match=None):
		"""Set qualities to None if there are no quality values"""
		self.name = name
		self.sequence = sequence
		self.qualities = qualities
		self.second_header = second_header
		self.match = match
		if qualities is not None and len(qualities) != len(sequence):
			rname = _shorten(name)
			raise FormatError("In read named {0!r}: length of quality sequence ({1}) and length "
				"of read ({2}) do not match".format(
					rname, len(qualities), len(sequence)))
	
	def __getitem__(self, key):
		"""slicing"""
		return self.__class__(
			self.name,
			self.sequence[key],
			self.qualities[key] if self.qualities is not None else None,
			self.second_header,
			self.match)

	def __repr__(self):
		qstr = ''
		if self.qualities is not None:
			qstr = ', qualities={0!r}'.format(_shorten(self.qualities))
		return '<Sequence(name={0!r}, sequence={1!r}{2})>'.format(_shorten(self.name), _shorten(self.sequence), qstr)

	def __len__(self):
		return len(self.sequence)

	def __richcmp__(self, other, int op):
		if 2 <= op <= 3:
			eq = self.name == other.name and \
				self.sequence == other.sequence and \
				self.qualities == other.qualities
			if op == 2:
				return eq
			else:
				return not eq
		else:
			raise NotImplementedError()

	def __reduce__(self):
		return (Sequence, (self.name, self.sequence, self.qualities, self.second_header,
		    self.match))


class FastqReader(SequenceReader):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.
		"""
		super(FastqReader, self).__init__(file)
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Yield Sequence objects
		"""
		cdef int i = 0
		cdef int strip
		cdef str line, name, qualities, sequence, name2
		sequence_class = self.sequence_class

		it = iter(self._file)
		line = next(it)
		if not (line and line[0] == '@'):
			raise FormatError("Line {0} in FASTQ file is expected to start with '@', but found {1!r}".format(i+1, line[:10]))
		strip = -2 if line.endswith('\r\n') else -1
		name = line[1:strip]

		i = 1
		for line in it:
			if i == 0:
				if not (line and line[0] == '@'):
					raise FormatError("Line {0} in FASTQ file is expected to start with '@', but found {1!r}".format(i+1, line[:10]))
				name = line[1:strip]
			elif i == 1:
				sequence = line[:strip]
			elif i == 2:
				if line == '+\n':  # check most common case first
					name2 = ''
				else:
					line = line[:strip]
					if not (line and line[0] == '+'):
						raise FormatError("Line {0} in FASTQ file is expected to start with '+', but found {1!r}".format(i+1, line[:10]))
					if len(line) > 1:
						if not line[1:] == name:
							raise FormatError(
								"At line {0}: Sequence descriptions in the FASTQ file don't match "
								"({1!r} != {2!r}).\n"
								"The second sequence description must be either empty "
								"or equal to the first description.".format(i+1,
									name, line[1:]))
						second_header = True
					else:
						second_header = False
			elif i == 3:
				if len(line) == len(sequence) - strip:
					qualities = line[:strip]
				else:
					qualities = line.rstrip('\r\n')
				yield sequence_class(name, sequence, qualities, second_header=second_header)
			i = (i + 1) % 4
		if i != 0:
			raise FormatError("FASTQ file ended prematurely")


class FastqReader2(SequenceReader):
	"""
	Reader for FASTQ files. Does not support multi-line FASTQ files.
	"""
	def __init__(self, file, sequence_class=Sequence):
		"""
		file is a filename or a file-like object.
		If file is a filename, then .gz files are supported.
		"""
		super(FastqReader2, self).__init__(file)
		self.sequence_class = sequence_class
		self.delivers_qualities = True

	def __iter__(self):
		"""
		Yield Sequence objects
		"""
		cdef bytearray buf = bytearray(1048576)
		cdef char[:] buf_view = buf
		cdef char* c
		cdef int endskip
		cdef bytes name_encoded
		cdef Py_ssize_t bufstart, bufend, pos, record_start, sequence_start
		cdef Py_ssize_t second_header_start, sequence_length, qualities_start
		cdef Py_ssize_t second_header_length, name_length
		cdef Py_ssize_t line
		readinto = self._file.buffer.readinto

		bufstart = 0
		line = 1
		# Read the input file in blocks
		while True:
			bufend = readinto(buf_view[bufstart:]) + bufstart
			if bufstart == bufend:
				# End of file
				break

			# Process this block. Parse one FASTQ record per iteration
			c = buf
			pos = 0

			# TODO
			# - an incomplete FASTQ is not seen as an error

			record_start = 0
			while True:
				if pos == bufend:
					break

				# Parse the name
				if c[pos] != '@':
					raise FormatError("Line {} in FASTQ file is expected to "
						"start with '@', but found {!r}".format(line, chr(c[pos])))
				pos += 1
				while pos < bufend and c[pos] != '\n':
					pos += 1
				if pos == bufend:
					break
				endskip = 1 if c[pos-1] == '\r' else 0
				name_length = pos - endskip - record_start - 1
				name_encoded = (<char*>buf)[record_start+1:pos-endskip]
				pos += 1
				line += 1

				# Parse the sequence
				sequence_start = pos
				while pos < bufend and c[pos] != '\n':
					pos += 1
				if pos == bufend:
					break
				endskip = 1 if c[pos-1] == '\r' else 0
				sequence = (<char*>buf)[sequence_start:pos-endskip].decode('ascii')
				sequence_length = pos - endskip - sequence_start
				pos += 1
				line += 1

				# Parse second header
				second_header_start = pos
				if pos == bufend:
					break
				if c[pos] != '+':
					raise FormatError("Line {} in FASTQ file is expected to "
						"start with '+', but found '{!r}'".format(line, chr(c[pos])))
				pos += 1
				while pos < bufend and c[pos] != '\n':
					pos += 1
				if pos == bufend:
					break
				line += 1
				endskip = 1 if c[pos-1] == '\r' else 0
				second_header_length = pos - endskip - second_header_start - 1
				if second_header_length == 0:
					second_header = False
				else:
					if (name_length != second_header_length or
							strncmp(<char*>buf+second_header_start+1,
								<char*>name_encoded, second_header_length) != 0):
						raise FormatError(
							"At line {}: Sequence descriptions in the "
							"FASTQ file don't match ('{}' != '{}').\n"
							"The second sequence description must be either "
							"empty or equal to the first description.".format(
								line, name_encoded.decode('ascii'),
								bytes(
									buf[second_header_start+1:pos-endskip]
								).decode('ascii')))
					second_header = True
				pos += 1
				line += 1

				# Parse qualities
				qualities_start = pos
				while pos < bufend and c[pos] != '\n':
					pos += 1
				if pos == bufend:
					break
				endskip = 1 if c[pos-1] == '\r' else 0
				qualities = (<char*>buf)[qualities_start:pos-endskip].decode('ascii')
				if pos - endskip - qualities_start != sequence_length:
					raise FormatError("At line {}: Length of sequence and "
						"qualities differ.".format(line))
				pos += 1
				line += 1

				yield Sequence(name_encoded.decode('ascii'), sequence, qualities)
				record_start = pos

			if pos == bufend:
				bufstart = bufend - record_start
				buf[0:bufstart] = buf[record_start:bufend]
