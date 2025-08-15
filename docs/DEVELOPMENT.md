## Development Environment Setup

### Required Tools

- PyCharm (recommended IDE)
- Python version and other package versions are written in 'environment.yml'

#### PyCharm Configuration

- Note: Do not set the directory as `Source Root`.
    - Reason: Setting it as `Source Root` causes import completions to start from one level below the package, making
      the test code unable to reference the intended modules properly.

---

## Contribution Policy

While this repository is not actively maintained, contributions are highly appreciated as they may help other users.
You are encouraged to:

- Submit questions or report bugs through the issue tracker.
- Propose fixes for reported issues via pull requests.

---

## Running Tests

```shell
python -m unittest discover -s tests -p "test_*.py"
```

