from unittest import TestCase

import ternpy


class TestJoke(TestCase):
    def test_is_string(self):
        s = ternpy.joke()
        self.assertTrue(isinstance(s, basestring))
