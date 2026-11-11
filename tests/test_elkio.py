"""Reader tests against files written here by hand in the documented binary layouts."""
import struct

import numpy as np

from exshg.elkio import fortran_records, read_evbse


def _record(payload, more=False):
    n = len(payload)
    head = -n if more else n
    return struct.pack("<i", head) + payload + struct.pack("<i", n)


def test_single_and_split_records(tmp_path):
    a = np.arange(5, dtype=np.int32).tobytes()
    b = np.linspace(0, 1, 7).tobytes()
    f = tmp_path / "x.bin"
    # second record split into two subrecords, as gfortran does above 2 GiB
    f.write_bytes(_record(a) + _record(b[:24], more=True) + _record(b[24:]))
    recs = fortran_records(f)
    assert len(recs) == 2
    assert np.array_equal(np.frombuffer(recs[0], np.int32), np.arange(5))
    assert np.allclose(np.frombuffer(recs[1], np.float64), np.linspace(0, 1, 7))


def test_read_evbse_fortran_order(tmp_path):
    nm = 3
    om = np.array([0.1, 0.2, 0.3])
    x = (np.arange(9) + 1j * np.arange(9)[::-1]).reshape(nm, nm)        # x[a, S]
    (tmp_path / "EVBSE.OUT").write_bytes(
        _record(struct.pack("<i", nm)) + _record(om.tobytes()) + _record(x.T.copy().tobytes()))
    om2, x2 = read_evbse(tmp_path)
    assert np.allclose(om2, om)
    assert np.allclose(x2, x)
