from flask import Flask, render_template

app = Flask(__name__)


@app.route('/', methods=['GET', 'POST'])
def hello_world():  # put application's code here
    return render_template('calculator.html', month="January")

# @app.route('/calculator')
# def calculator():


if __name__ == '__main__':
    app.run()
